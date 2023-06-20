function [ x ] = IterGSS_improve(y, rho, K, Method, Nit, init)
% This function realizes iterative generalized structured shrinkage based
% on TQWT
%% Input %%%%%%%%%%
%   y      : noisy data
%   A, AH  : function handles for A and its transpose
%   normA  : the norm of the transformation
%   rho    : rho >= maximum eigenvalue of A'A
%   K      : K-sparsity parameters
%   Method : a structure with information of the method
%   Nit    : number of iterations
%   init   : the inital point of x
%% Output %%%%%%%%%%
%   x      : Recovery

% Reference: 'Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware',
% IEEE Transactions on Instrumentation and Measurement, 2020
% https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2019.6

% initialization
AHy = (y);
if ~exist('init', 'var')
    init = AHy;
end
if ~exist('Nit', 'var')
    Nit = 50;
end
% cost = zeros(Nit , 1);
mu = 1 / rho;
x = init;
x_old = init;
Ax = x;
iter = 1;

while iter <= Nit
%     fprintf(['Iteration:' num2str(iter) '\n'])
    % forward operator
    % x = x - mu * (A^T(A(x) - y))
    tmp = x; 
    AHAx = (Ax); 
%     for i = 1:numel(tmp) 
        tmp = tmp + mu * ( AHy - AHAx);  
%     end
     
    % backward (thresholding)
    % x = soft(x, mu * lam);
    Temp = [];
%     for i = 1:numel(x)
        Temp = [Temp ; tmp(:)];
%     end
    c = K_sparsity(Temp, K);
    for i = 1 %numel(x):-1:1
        if ~strcmp(Method.Name, 'WGL')
            Method.lambda = c ;
            Phi = Shrinkage(Method.Name, Method);
            x = Phi(tmp);               
            
        else
            if i == 1 %numel(x)
                Size = Method.Initial_Size;
                x = generalized_structured_shrinkage(tmp, c , Size ,  Method);
            else
                Length_Before = length(tmp{i+1});
                Length_Now = length(tmp{i});
                if Length_Before ~= Length_Now
                    Size = (Size - 1) * 2 + 1;
                end
                x{i} = generalized_structured_shrinkage(tmp{i}, c , Size, Method);
            end
        end
    end    
 
    % reconstruction
    x=x(:);
    Ax = x;
    
    stop = 0;
    total = 0;
     for i = 1 %numel(x)
        stop = stop + sum((x - x_old).^2);
        total = total + sum((x).^2);
    end   

    if sqrt(stop) / sqrt(total) < 1e-6
        break;
    end
    iter = iter + 1;
    x_old = x;
end

111;

