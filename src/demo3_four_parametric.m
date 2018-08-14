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
% - demonstration of four-parametric sharp Loewner majorants
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

n_points    = 2000;     % number of points for the discrete value for B_2
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

%% split the matrix
A2  =   A( 1:n, 1:n );
a   =   A( 1:n, n+1 );
a0  =   A( n+1, n+1 );

%% majorant for the upper d x d - block matrix
[ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A2, p0, tol );
% if case 2, then mu_hat is not contained in the admissible domain
if( idx_case == 2 )
    fprintf('old mu_hat: %15.6e\n', mu_hat);
    mu_hat = mu_hat + 1e-5*(max(Lambda)-min(Lambda));
    fprintf('new mu_hat: %15.6e\n', mu_hat);
end

% ATTENTION: cutoff must be a high number, otherwise for beta_1 --> mu_hat
%            the beta_2 value can be too small (since limited to b2cutoff)
[ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, 2*mu_hat-Lambda(end), Lambda, p, idx_case, n_points, b2cutoff );

MajorantInfoP2( Lambda, p, idx_case );
PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case );

beta3vec = zeros(size(beta1vec));
beta4vec = zeros(size(beta1vec));

D   = @(b1) b1*eye(n) - diag(Lambda);
l   = U'*a;
for i=1:length(beta2vec)
    b1 = beta1vec(i);
    b2 = beta2vec(i);

    %% now check for the possibility of gathering further special cases:
    % if Delta A is not singular, then for any given zeta_0
    % k0 can be provided in closed form:
    Di = pinv(D(b1));
    nn  = Di*p;
    Ai = pinv(D(b1)+b2*p*p');

    % optimal zeta_0:
    b4      = ( l' * nn ) / ( p'* nn );
    z       = b4*p - l;
    b3      = a0 + z' * (Ai*z);

    beta3vec(i) = b3;
    beta4vec(i) = b4;

%     %% OPTIONAL: investigate the first few eigenvalues
%     % of the difference matrix:
%     B4  = [ b1*eye(n) + b2*p0*p0', b4 * p0; b4 * p0', b3 ];
%     D4  = B4 - A; %[ diag(Lambda), zeta_tilde; zeta_tilde' , k ];
%     e   = eig( D4 );
%     e   = sort(e,'ascend');
%     fprintf('    smallest three eig. of B4 - A   %14.6e  %14.6e  %14.6e\n', e(1), e(2), e(3) );
end

DATA = [    beta1vec; ...
            beta2vec; ...
            beta3vec; ...
            beta4vec ];
dlmwrite('b4_example.dat', DATA', 'precision', '%16.8e', 'delimiter', ' ');


%% Semidefinite programming problem
% TODO:
% objective function
% 