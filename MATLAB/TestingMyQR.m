clear; close all; clc;

% Testing my QR pivoting algo
A = rand(6,8);
[Q,R,P] = qr_householderpivotingIII(A)
[~,~,P] = qr(A, 'vector')
%error = norm(Q*R - A, 2)