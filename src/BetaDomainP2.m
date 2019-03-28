function [ beta1vec, beta2vec ] = BetaDomainP2( b1_min, b1_max, Lambda, p, idx_case, n_points, b2cutoff )
%% BETADOMAINP2 computes the boundary of the critical beta1/2 domain
% for given outputs of the (previously called) SHARPLOEWNERMAJORANTP2.
% The outputs can be plotted using PLOTCRITICALB2 or they can be processed
% in discrete optimization quite easily.
%
% For subsequent processing for set majorization use INTERSECTMAJORANTS.
%
% INPUTS
% b1_min        upper bound for beta_1 ([] --> mu_hat )
% b1_max        upper bound for beta_1 ([] --> mu_hat + 3*(Lmax-Lmin)
% Lambda        d-vector containing the (descendingly sorted) eigenvalues
%               of the symmetric matrix A under consideration
% p             normalized d-vector in the coordinate system of A
% idx_case      case to be considered (1 --> p = e_i , 2 --> p inclined)
% n_points      [OPTIONAL] number of outputs (default: 100)
% b2cutoff      [OPTIONAL] cutoff for beta_2 (default: 4 * mu_hat - mu_p)
%
% OUTPUTS
% beta1vec      row-vector of length n containing beta_1 values
% beta2vec      row-vector of length n containing beta_2 values
%
% See also SHARPLOEWNERMAJORANTP2, PLOTCRITICALB2, INTERSECTMAJORANTS.
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
% https://github.com/EMMA-Group/LoewnerMajorant
%
% CHANGELOG
% Aug 08, 2018      initial creation


tol     = 1e-8;
m       = length(p);
mu_p    = p'*(Lambda.*p);
Q       = orth(eye(m)-p*p');
AQ      = Q'*diag(Lambda)*Q;
mu_hat  = max(eig(AQ));

%% set b1_min = mu_hat if nothing is specified
if( ~exist('b1_min','var') )
    b1_min = mu_hat;
end

%% if b1_max is empty, then compute some default value that should be
% interesting, i.e. not too close to 'saturation'/linear b1-b2-relation
if( ~exist('b1_max','var') )
    b1_max  = b1_min + 3*(max(Lambda)-min(Lambda));
end
if( isempty(b1_max) )
    b1_max  = b1_min + 3*(max(Lambda)-min(Lambda));
end

%% test for plausibility of arguments
if( ~exist('n_points', 'var') )
    n_points = 100;
    fprintf('%% using default number of steps n=%i\n', n_points);
end

if( ~exist('b2cutoff', 'var') )
    b2cutoff = 4*mu_hat-mu_p;
end

if( mu_hat - b1_min > tol*abs(mu_hat))
    warning('WARNING: b1_min is too small (i.e. smaller than mu_hat)');
end

beta1vec = linspace(b1_min, b1_max, n_points);

%% CASE 1: p is coaxial to some eigenvectors of A
if( idx_case == 1 )
    % p0 is paralel to the eigenspace of A
    beta2vec = mu_p - beta1vec;
    return;
end

%% CASE 2: p is not coaxial to some eigenvectors of A

% short form of p' * (b1*I - Lambda)^-1 * p
factor      = @(b1) -1/sum( (p.^2)./(b1-Lambda) );

% evaluate beta_2 at *all* points
% (even the ones which need special care; see below below)
% calculation cf. Case 2.2
beta2vec    = arrayfun( factor, beta1vec );

% if p orthogonal to the eigenspace of Lambda(1) ???
J           = find( abs( Lambda - Lambda(1) ) < tol * norm(Lambda) );
orthogonal  = ( norm(p(J)) < tol );

% get the indices of the beta_1 values that match Lambda(1) [if any]
% for special treatment
I           = find( abs( beta1vec - Lambda(1) ) < tol * norm(Lambda) );
if( ~isempty(I))
    if( ~orthogonal )
        %% Case 2.1a
        fprintf('case 2.1a\n');
        beta2vec( I ) = 0.;
        if( abs( beta1vec(1) - mu_hat ) < tol )
            beta2vec(1) = b2cutoff;
        end
    else
        %% Case 2.1b: need the pseudo-inverse in select points
        for j=I
            beta2vec(j) = - 1/( p' * pinv( beta1vec(j) * eye(length(Lambda)) - diag(Lambda) ) * p );
        end
        fprintf('case 2.1b\n');
    end
else
    if( abs( beta1vec(1) - mu_hat ) < tol )
        beta2vec(1) = b2cutoff;
    end

end
beta2vec = min(beta2vec, b2cutoff);

end
