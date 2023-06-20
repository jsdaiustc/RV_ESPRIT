function Q=real_trans(N)

% J=zeros(floor(N/2));
% for ii=1:floor(N/2)
%    J(floor(N/2)+1-ii, ii)=1;
% end

J=flip(  speye(floor(N/2))  );
K=speye(floor(N/2));
X=sparse( N-2*floor(N/2), floor(N/2));
Y=sparse(floor(N/2), N-2*floor(N/2));
Q=1/sqrt(2)*[ K,            Y,                            J 
              X,   sqrt(2)*ones(1,N-2*floor(N/2)),        X
             1i*J,          Y,                          -1i* K ];
