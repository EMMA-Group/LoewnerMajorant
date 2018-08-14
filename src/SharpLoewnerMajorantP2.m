function [ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A, p0, tol )
%% compute the domain of admissible parameters beta1, beta2 that produce
% sharp Loewner majorants of a matrix A
%
% INPUTS
% A         symmetric n x n matrix
% p0        n dim. vector (in the standard basis)
% tol       [OPTIONAL]  a tolerance for the detection of colinearity of p and
%                       eigenvectors of A
%
% OUTPUTS
% U         orthonormal eigenbasis of A
% Lambda    eigenvalues of A (row vector)
% mu_p      mu_p = p.A.p
% mu_hat    max. eigenvalue of A in the vector space orthogonal to p
% idx_case  number of the case (can be used to automatically find
%           appropriate beta_2 leading to sharp bounds or to draw
%           the domain boundary)
%
% See also BETADOMAINP2, PLOTCRITICALB2
%
%
% related to the article
% Sharp Loewner majorants for a set of symmetric matrices
% by Felix Fritzen and Mauricio Fernandez
% [ submitted to AMCS (Aug 2018) ]
%
% CONTENT OF THIS FILE
% - examples for the two-parametric matrix majorization (se Fig. X, Y)
% - illustrates how to access the different tools
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
% XXXXXXXXXXXXXXXXXXXX
%
% CHANGELOG
% Aug 08, 2018      initial creation


%% set default tolerance, if needed
if( ~exist('tol', 'var') )
    tol = 1e-7;
    fprintf('%% using default tolerance; tol=%12.6e\n', tol);
end

%% check if tolerance is admissible, i.e. positive
if( tol <= 0 )
    error(sprintf('ERROR: tolerance must be nonnegative, but tol=%12.5e', tol));
end

%% check if A is admissible
m   = size(A,1);
if( size(A,2) ~= m )
    error(sprintf('ERROR: A must be symmetric and square, but received a %i x %i matrix', m, size(A,2)));
end
if( norm( A-A', 'fro') > tol*norm(A,'fro') )
    error('ERROR: A must be symmetric, but the skew part has Frobenius norm larger than TOL * ||A||_FRO');
end

%% check if p0 is admissible
if( size(p0,1) ~= m)
    error('ERROR: dimension of p0 must match the one of A');
end
if( size(p0,2) ~= 1)
    error('ERROR: p0 is exepected as a column vector (i.e. only one single column)');
end
if( abs(1-norm(p0)) > tol )
    error('ERROR: p0 must be normalized');
end


%% setup phase:
%  compute eigenbasis
%  perform coordinate transform
%  decide on which case is occurring
[ U, L ]    = eig(A);
Lambda      = diag(L);
[Lambda, I] = sort(Lambda, 'descend');
U           = U(:,I);
L           = diag(Lambda);
p           = U' * p0;
Q           = orth(eye(m)-p*p');
AQ          = Q'*L*Q;

%% check if p is co-linear to some eigenvector
mu_p        = p'*L*p;
mu_hat      = max(eig(AQ));

%% figure out if the eigensystem is coaxial to p0 (CASE 1)
[I]     = find( (abs(mu_p - Lambda)) < tol*sum(abs(Lambda)) );
if( length(I) >= 1 )
    %% p0 belongs to an eigenspace of A
    % fprintf('%% found p''*A*p to match an eigenspace of dimension %i\n',length(I));
    idx_case    = 1;
    return;
end

%% p is not parallel to an eigensystem (CASE 2)
idx_case    = 2;
end