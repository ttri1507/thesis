clc
clear

% example here: https://www.youtube.com/watch?v=-3my1TkyFiM
W = [217 125 88 109];
A = [58 44 26 23;
    25 29 12 17;
    43 25 23 29];
B = [120; 80; 95];
tic
cvx_solver mosek
cvx_begin
    variable X(4,1) binary;
    obj = W*X;
    maximize(obj)
    subject to
        A*X <= B;
cvx_end
toc
% % Input data
% m = 16; n = 8;
% As = randn(m,n);
% bs = randn(m,1);
% 
% % Matlab version
% x_ls = As \ bs;
% 
% % cvx version
% cvx_begin
%     variable x(n)
%     minimize( norm(As*x-bs) )
% cvx_end
