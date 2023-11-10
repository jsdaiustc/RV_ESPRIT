clear;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

range=[0,1000];      %%%%%%滤波范围 
M=1000;               %%%%%%窗口大小
fw=800;             %%%%%%%坐标横轴大小

%% Paderborn University dataset Fs=64000 理论故障频率123.3 Hz
 y=load('N15_M07_F04_KI18_10.mat');
 y=y.N15_M07_F04_KI18_10.Y(7).Data; 
 y=y';
 y=y(1:16384);
Fs=64000;
N=length(y);
y_env=hilbert(abs(hilbert(y))-mean(abs(hilbert(y))));


%% FFT
y_fft= fft(abs(hilbert(y))-mean(abs(hilbert(y))));
subplot(3,2,1);
plot([0:1:length(y)-1]*Fs/length(y), abs(y_fft)/(N/2)); 
axis([0 fw 0 0.2]);
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Times[s])');
ylabel('\fontname{Times New Roman}Amplitude\fontname{Times New Roman}(m/s^2)');
title('(a) FFT','fontname','Times New Roman');


%% GSL
F=([1:length(y)]-1)*Fs/length(y);
mu=GSL(y);
y_envo= fft(abs(hilbert(mu))-mean(abs(hilbert(mu))));
subplot(3,2,2);
plot(F, abs(y_envo)/(N/2));
axis([0 fw 0 0.2])
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
ylabel('\fontname{Times New Roman}Amplitude\fontname{Times New Roman}(m/s^2)');
title('(b) GSL');


%% BPD
rho = 1;
K2 = round(N*0.04);
Method2.Name = 'WGL';
Method2.Initial_Size = 5;
Method2.SubName = 'MC';
Method2.gamma = 2;
Method2.window = 'gausswin';
Z2= IterGSS_improve(y, rho, K2, Method2);
y_BPD= real(Z2);
y_AdaESPGL_BPD=abs(fft(abs(hilbert(y_BPD)) -mean(abs(hilbert(y_BPD)))  ))/(N/2);
subplot(3,2,3);
plot(F  ,y_AdaESPGL_BPD)
axis([0 fw 0 0.2])
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
ylabel('\fontname{Times New Roman}Amplitude\fontname{Times New Roman}(m/s^2)');
title('(c) BPD')



%% GSSA
Q = 2;
r = 5;
J =10;
now = ComputeNow(N,Q,r,J,'radix2');
AH = @(Sig) tqwt_radix2(Sig, Q, r, J);
A = @(w) itqwt_radix2(w, Q, r , N);
lam = 1.0 * now;
rho = 1;
load Performance_Comparison_Combination_K_Index_Size5_Sigma6.mat
K1 = round(N*0.04); 
Method1.Name = 'WGL';
Method1.Initial_Size = 5;
Method1.SubName = 'MC';
Method1.gamma = 2;
Method1.window = 'gausswin';
z1 = IterGSS(y, A, AH, lam, rho, K1, Method1);
P_GSSA = real(A(z1));
y_GSSA_enve=abs(fft(abs(hilbert(P_GSSA)) -mean(abs(hilbert(P_GSSA)))  ))/(N/2);
subplot(3,2,4);
plot(F  ,y_GSSA_enve)
axis([0 fw 0 0.2])
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
ylabel('\fontname{Times New Roman}Amplitude\fontname{Times New Roman}(m/s^2)');
title('(d) GSSA')


%% TLS_ESPRIT
[F_est_all,~,s_p] =TLS_ESPRIT(Fs,y_env,M);
figure(1);
subplot(3,2,5);
stem(F_est_all, abs(s_p),'marker','none');
axis([0 fw 0 0.2]);
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
ylabel('\fontname{Times New Roman}Amplitude\fontname{Times New Roman}(m/s^2)');
title('(e) TLS-ESPRIT','fontname','Times New Roman');



%% RV_ESPRIT
[hat_f,~,hat_s]=RV_ESPRIT(y_env,M,Fs,range);
figure(1);
subplot(3,2,6);
stem(hat_f,abs(hat_s),'marker','none');
axis([0 fw 0 0.2]);
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
ylabel('\fontname{Times New Roman}Amplitude\fontname{Times New Roman}(m/s^2)');
title('(f) RV-ESPRIT','fontname','Times New Roman');
