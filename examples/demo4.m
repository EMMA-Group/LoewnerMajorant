clc;
clear all;
close all;
addpath(genpath('..'))

%% Nonlinear function, matrix A4 \in S_{n+1} and p0 \in R^{n}

% Nonlinear function to be minimized: Frobenius norm of four-parametric form
phi = @(beta,n) sqrt(n*beta(1)^2 + beta(2)^2 + beta(3)^2 + 2*beta(4)^2);

% Load A4 and p0
load('data/A4_p0.mat');
n = length(p0);
disp('Eigenvalues of A4');
eig(A4)

% Information on case for A2 part of A4 matrix
B4Preparations(A4,p0,1);

% B4 in original space and eigenvalues of difference matrix
B40 = @(b) [b(1)*eye(n) + b(2)*p0*p0',b(4)*p0; b(4)*p0',b(3)];
eigC4 = @(b) sort(eig(B40(b) - A4));

%% Compute minimum along critical b1-path

disp('**********************************************************************')
disp('Compute minimum along critical b1-path')

[bmin1,phimin1] = B4minSharpPath(@(b) phi(b,n),A4,p0)

disp('Lowest 4 eigenvalues of difference matrix C4')
ev = eigC4(bmin1)';
ev(1:4)

%% Numerical optimization

disp('**********************************************************************')
disp('Numerical optimization and eigenvalues of difference matrix C4')

% Preparations and shifted mu avoiding case differentiation for dim(ker(C4))>1
[~,Lv,l,l0,p,~,mu_hat] = B4Preparations(A4,p0);
if mu_hat > 0
    mu_shifted = 1.01*mu_hat;
else
    mu_shifted = 0.99*mu_hat;
end

% Nonlinear condition function based on lemma 4 and lemma 5 of manuscript
cond = @(x) B4ShiftedConditions(x,Lv,l,l0,p);

% Initial guess in the vicinity of bmin1 slightly toward trivial majorant
Lmax4 = max(eig(A4));
trivial = [Lmax4,0,Lmax4,0];
bmin20 = bmin1+0.01*(trivial-bmin1);

% Optimize numerically with mu_shifted<=b1 and shifted conditions
[bmin2,phimin2] = fmincon(@(x) phi(x,n),bmin20,[],[],[],[],[mu_shifted;-inf;-inf;-inf],[],cond)

disp('Lowest 4 eigenvalues of difference matrix C4')
ev = eigC4(bmin2)';
ev(1:4)