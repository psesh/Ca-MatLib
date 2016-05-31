clear; close all; clc;

% Testing my QR pivoting algo
A = rand(6,4);
[Q,R,P] = qr_householderpivotingIII(A)

%error = norm(Q*R - A, 2)