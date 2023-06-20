function [x_f,s_p,hat_s]=RV_ESPRIT(y,M,Fs,range)
N=length(y);
%%
f=fft(y);
fs=((1:N)-1)*Fs/N;
f(fs>range(2))=0;
y=ifft(f);

%% spatial smoothing
G=N+1-M;
x=y(1:M);
for kk=1:G-1
    x=[x,  y(kk+1:kk+M)] ;
end


Q_M=real_trans(M);
x_trans=Q_M*x;
x_real=real(x_trans);

Y=x_real*x_real';
[V,T]=eig( Y );
lam=diag(T);
[~,I]=sort(lam,'descend');%特征降序排列
V=V(:,I);%特征向量排序

Q_M1=real_trans(M-1);
H1=2*real(Q_M1*Q_M(:,1:end-1)');
H2=2*imag(Q_M1*Q_M(:,2:end)');

L=10;
Us=V(:,1:L);
S1=H1*Us;
S2=H2*Us;
S=[S1,S2];
[P,~]=svd(S'*S);
P12=P(1:L,L+1:2*L);
P22=P(L+1:2*L,L+1:2*L);
Phi=-P12*inv(P22);
mu=real(eig(Phi));
RS_angle=2*atan(mu);
IND_negative=find(RS_angle<0);
RS_angle(IND_negative)=RS_angle(IND_negative)+ 2*pi; 
hat_f=RS_angle/(2*pi/Fs);
hat_A=exp(2*pi*1j*hat_f*(0:N-1)/Fs).';
hat_s=pinv(hat_A'*hat_A)*hat_A'*y;
s_p=mean(abs(hat_s).^2,2);
ind_remove=find(s_p< mean(s_p)/20);
s_p(ind_remove)=[];
hat_f(ind_remove)=[];

%% 去掉故障频率中的倍频%%%%
 [fn,~,F_active,~]=search_Pm(s_p,hat_f); 
 hat_f=F_active;
 hat_f=[fn;hat_f];
 
  
 %% 求模差，并选取差值最小的故障频率%%%%
for ii=1:length(hat_f)
    F_new= (hat_f(ii)*[1:1:10])';
    hat_A=exp(2*pi*1j*F_new*(0:N-1)/Fs).';
    hat_s=inv(hat_A'*hat_A)*hat_A'*y;
    s_p=mean(abs(hat_s).^2,2); 
    hat_A=exp(2*pi*1j*F_new*(0:N-1)/Fs).';
    hat_s=inv(hat_A'*hat_A)*hat_A'*y;
    e=y-hat_A*hat_s;
    disp(ii)=e'*e/N;
    norm_e(ii)=norm(e)^2/norm(y)^2;
   
    
    s_p=sqrt(mean(abs(hat_s).^2,2));        
    xx{ii}=F_new;
    yy{ii}=s_p;    
end


ind_opt= 1;
[~,F_location]=min(disp);
if disp(ind_opt)>   disp(F_location)*1.01
    ind_opt=F_location;
end
x_f=xx{ind_opt};

if norm_e(F_location)>0.92
    s_p=s_p*1e-10;
else
    s_p=yy{ind_opt};
end

hat_s=s_p;



