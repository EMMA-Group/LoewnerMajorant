clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% related to the article
% Sharp Loewner majorants for a set of symmetric matrices
% by Felix Fritzen and Mauricio Fernandez
% [ submitted to AMCS (Aug 2018) ]
%
% CONTENT OF THIS FILE
% - demonstration of four-parametric optimization with
%   nonlinear contraints for set of matrices
%
% LICENSE
% This software is distributed under the GNU GPL v3.
% An active contribution to this software is welcome (see CONTACT below).
%
% Further licensing information is contained in the file 'license'
% distributed with this code.
%
% CONTACT
% Felix Fritzen        felix.fritzen@mechbau.uni-stuttgart.de
% Mauricio Fernandez   mauricio.fernandez@mechbau.uni-stuttgart.de
% http://www.mechbau.uni-stuttgart.de/EMMA
%
% GITHUB
% This software is hosted on GITHUB. It is accessible via
% https://www.github.com/EMMA-Group/LoewnerMajorant
%
% CHANGELOG
% Aug 14, 2018      initial creation
%

n_points    = 4000;     % number of points for the discrete value for B_2
n           = 10;       % dimension of the matrix
tol         = 1e-8;
b2cutoff    = 1e30;     % very high number since otherwise,
%                         beta1-->mu_hat leads to numerical
%                         issues in the validation

%% generate random matrix
A   = randn(n+1,n+1);
A   = A + A' - 0.25*n*eye(n+1); % make the matrix possibly indefinite
%% define random p0
p0  = randn(n,1);
p0  = p0 / norm(p0);

%% store data:
% dlmwrite('../data/paper_b4_matrix_A.dat',A,'precision','%16.8e','delimiter',' ');
% dlmwrite('../data/paper_b4_vector_p0.dat',p0,'precision','%16.8e','delimiter',' ');

%% use data from the article:
 A  = dlmread('../data/paper_b4_matrix_A.dat');
 p0 = dlmread('../data/paper_b4_vector_p0.dat');

%% split the matrix
A2  =   A( 1:n, 1:n );
a   =   A( 1:n, n+1 );
a0  =   A( n+1, n+1 );

%% majorant for the upper n x n - block matrix
[ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A2, p0, tol );

% transform a into l
l   = U'*a;

% if case 2, then mu_hat is not contained in the admissible domain
if( idx_case == 2 )
    fprintf('%% old mu_hat: %15.6e\n', mu_hat);
    mu_hat = mu_hat + 1e-5*(max(Lambda)-min(Lambda));
    fprintf('%% new mu_hat: %15.6e\n', mu_hat);
end

%% evluate critical domain of the upper n x n block of the difference matrix
% * ATTENTION: cutoff must be a high number, otherwise for beta_1 --> mu_hat
% *            the beta_2 value can be too small (since limited to b2cutoff)
[ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, 2*mu_hat-Lambda(end), Lambda, p, idx_case, n_points, b2cutoff );

%% optional output
% MajorantInfoP2( Lambda, p, idx_case );
% PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case );

beta3vec = zeros(size(beta1vec));
beta4vec = zeros(size(beta1vec));

%% compute double-sharp Loewner Majorants:
% set beta1, beta2 to critical values and optimize
% beta4, beta3 to get the smallest possible beta3
% a double sharp majorant yielding.
% (for validation, check the smallest eigenvalues of the difference matrix
% D4 = [ D(b1) + b2*p*p', b4*p - l; b4*p'-l; b3-a0 ];
% which should be zero up to numberical precision)

D   = @(b1) diag(b1-Lambda);
for i=1:length(beta2vec)
    b1 = beta1vec(i);
    b2 = beta2vec(i);

    %% now check for the possibility of gathering further special cases:
    % if Delta A is not singular, then for any given zeta_0
    % k0 can be provided in closed form:
    Di  = pinv(D(b1));
    nn  = Di*p;
    Ai  = pinv(D(b1)+b2*p*p');

    % optimal zeta_0:
    b4      = ( l' * nn ) / ( p'* nn );
    c       = b4*p - l;
    b3      = a0 + c' * (Ai*c);

    beta3vec(i) = b3;
    beta4vec(i) = b4;

end
%% gather all beta_i in one 4 x n_points matrix
DATA = [    beta1vec; ...
            beta2vec; ...
            beta3vec; ...
            beta4vec ];
%  dlmwrite('../data/paper_B4_example_beta1234.dat', DATA', 'precision', '%16.8e', 'delimiter', ' ');


%% Semidefinite programming problem

% objective function: Frobenius norm squared
phi = @(beta) sqrt(n*beta(1)^2 + beta(2)^2 + beta(3)^2 + 2*beta(4)^2);

% compute discrete values of phi along the double sharp majorants
% computed earalier
F = zeros(1, length(beta1vec));
for i=1:length(beta1vec)
    F(i) = phi(DATA(:,i));
end

% plot phi over beta1
plot(beta1vec, F, '-', 'linewidth', 2, 'color', 'black');
axis([mu_hat-1 20 10 50]);

% get optimum
[Fmin, I]   = min(F);
bmin        = DATA(:,I);
bmin0       = bmin;
Fmin0       = Fmin;

%% use constraint optimization to find (possibly sharp) majorants minimizing
% the objective function using numerical optimization

%% define some options for FMINCON
options = optimoptions('fmincon');

% Set OptimalityTolerance to 1e-3
options = optimoptions(options, 'OptimalityTolerance', 1e-12); 

% Set the Display option to 'iter', the StepTolerance
% and the constraint tolerance (!!! important !!!)
options.Display = 'iter';
options.StepTolerance = 1e-12;
options.ConstraintTolerance = 1e-12;


%% define constraints (here: pass needed parameters to the function)
C4      = @(x) [ diag(x(1)-Lambda)+x(2)*p*p', x(4)*p-l; x(4)*p'-l', x(3)-a0 ];
cond    = @(x) conditions( x,Lambda, p, l, a0);
phi     = @(x) Phi_Frobenius(x, n);
bmin    = fmincon( phi, bmin0 + 0.05*[rand(2,1);zeros(2,1)], [], [], [], [], [mu_hat;-inf;-inf;-inf], [], cond, options );
Fmin    = phi(bmin);
fprintf('\n\n**********************************************************************************************\n');
fprintf('* optimal value for double sharp Loewner Majorants   %16.8e\n', Fmin0);
fprintf('*          beta = [ %12.6e; %12.6e; %12.6e; %12.6e ];\n', ...
    bmin0(1), bmin0(2), bmin0(3), bmin0(4) );
e = sort( eig(C4(bmin0)), 'ascend');
fprintf('* two smallest eigenvalues of the difference matrix: %12.5e, %12.5e\n', e(1), e(2));
fprintf('**********************************************************************************************\n\n\n');
fprintf('**********************************************************************************************\n');
fprintf('* optimal value for all Loewner Majorants            %16.8e\n', Fmin);
fprintf('*          beta = [ %12.6e; %12.6e; %12.6e; %12.6e ];\n', ...
    bmin(1), bmin(2), bmin(3), bmin(4) );

e = sort( eig(C4(bmin)), 'ascend');
fprintf('* two smallest eigenvalues of the difference matrix: %12.5e, %12.5e\n', e(1), e(2));
fprintf('**********************************************************************************************\n\n\n');
fprintf('**********************************************************************************************\n');
fprintf('* change minimum over double sharp majorants         %16.8e (relative: %10.5f)\n', Fmin-Fmin0, Fmin/Fmin0-1);
dbmin=bmin-bmin0;
fprintf('* change of beta  [ %12.6e, %12.6e, %12.6e, %12.6e ]\n', dbmin(1), dbmin(2), dbmin(3), dbmin(4) );
fprintf('**********************************************************************************************\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return Forbenius norm of K0
function [ phi ] = Phi_Frobenius(x, n)
phi     = sqrt( n*x(1)^2+x(2)^2+x(3)^2+2*x(4)^2 );

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq]=conditions(x, Lambda, p, l, a0)
%% nonlinear inequality constraints:
% beta2 >= beta2_c(beta1)
% beta3 >= a0 + (beta3*p-l)' * (B2(beta1,beta2)-Lambda)^-1 * (beta3*p-l)
% and check if beta4 is admissible in the case of critical beta1/2-pairs
tol = 1e-14;
ceq = [];

n   = p./(x(1)-Lambda);
c   = [ - x(2) - 1 / (n'*p); 1 ];
z   = x(4)*p - l;

%% check if beta2=beta2_crit(beta1) --> C has a nullspace!
if( c(1)>-tol )
    %% check if beta4 is admissible (i.e. such that a Loewner Majorant is
    % still obtained)
    b4crit  = n'*l/(n'*p);
    d       = abs(b4crit-x(4));
    if( d > tol )
        %% beta4 is inadmissible --> penalty function
        fac  = 1e12;
        %% Attention: fac and 1/tol should not differ by more than
        % sqrt(accuracy_fmincon)
        c(2) = d * fac;
        return;
    else
        %% check if beta3 is ok
        zi  = pinv( diag(x(1)-Lambda) + x(2)*p*p') * z;
        c(2)= a0 + z'*zi - x(3);
        return;
    end
end

%% check if constraint for beta3 is fulfilled
% (note: zi computes exactly since the matrix is positive definite)
z   = x(4)*p - l;
zi  = ( diag(x(1)-Lambda) + x(2)*p*p') \ z;
c(2)= a0 + z'*zi - x(3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
