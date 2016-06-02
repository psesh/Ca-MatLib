clear; close all; clc;

% Testing my QR pivoting algo
A = rand(6,8);
[Q,R,P] = qr_mgs_pivoting(A)A *
[~,~,P] = qr(A, 'vector')