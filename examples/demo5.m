clc;
clear all;
close all;
addpath(genpath('..'))

%% Nonlinear function, matrix A4 \in S_{n+1} and p0 \in R^{n}

% Nonlinear function to be minimized: Frobenius norm of four-parametric
% form
phi = @(beta,n) sqrt(n*beta(1)^2 + beta(2)^2 + beta(3)^2 + 2*beta(4)^2);

% Load matrix set
load('data/A4set_p0.mat','A4','p0');
N = length(A4);
n = length(p0);
fprintf('Number of matrices: %i\nDimension of p0: %i\n',N,n);

% B4 in original space and eigenvalues of difference matrix
B40 = @(b) [b(1)*eye(n) + b(2)*p0*p0',b(4)*p0; b(4)*p0',b(3)];
eigC4 = @(b,A4) sort(eig(B40(b) - A4));

% Get cases for each matrix, just out of interest
disp('Case for each matrix:')
B4Preparations(A4,p0,1);

%% Minimization along sharp b1-paths with b4=0

disp('**********************************************************************')
disp('Compute minimum along sharp b1-paths with b4=0')

[bmin1,phimin1] = B4minSharpPath(@(b) phi(b,n),A4,p0)

% Check eigenvalues of difference matrices
disp('Check 4 lowest eigenvalues of difference matrices')
for i=1:length(A4)
    ev = eigC4(bmin1,A4{i})';
    disp(i)
    disp(ev(1:4))
end

%% Numerical optimization

disp('**********************************************************************')
disp('Numerical optimization')

% Preparations and shifted mu avoiding case differentiation for 
% dim(ker(C4))>1
[~,Lv,l,l0,p,mu_P,mu_Q,c_Q] = B4Preparations(A4,p0);
mu_Q_max = max(cell2mat(mu_Q));
mu_Q_s = mu_Q_max+0.01*abs(mu_Q_max);

% Nonlinear condition function based on lemma 4 and lemma 5 of manuscript
% for mu_Q_s ot be used later on (i.e., mu_Q_max < mu_Q_s <= b1)
cond = @(x) B4ShiftedConditions(x,Lv,l,l0,p);

% Initial guess in the vicinity of bmin1 slightly towards trivial majorant
Lmax4 = max(arrayfun(@(i) max(eig(A4{i})),1:N));
trivial = [Lmax4,0,Lmax4,0];
bmin20 = bmin1+0.01*(trivial-bmin1);

% Optimize numerically with mu_shifted<=b1 and shifted conditions
[bmin2,phimin2] = fmincon(@(x) phi(x,n),bmin20,[],[],[],[],[mu_Q_s;-inf;-inf;-inf],[],cond)

% Check eigenvalues of difference matrices
disp('Check 4 lowest eigenvalues of difference matrices')
for i=1:length(A4)
    ev = eigC4(bmin2,A4{i})';
    disp(i)
    disp(ev(1:4))
end