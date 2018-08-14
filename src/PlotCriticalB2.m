function h = PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case, c )
%% PLOTCRITICALB2 displays the results as well as the asymptotic curves
% using the output of the routines SHARPLOEWNERMAJORANTP2 and BETADOMAINP2.
%
% INPUTS
% beta1vec      row-vector containing the beta_1 values to be plotted
% beta2vec      row-vector containing the beta_2 values to be plotted
% Lambda        diagonal representation of A
% p             p0 in the orthogonal eigensystem of A
% idx_case      case determined by SharpLoewnerMajorantP2
% c             [ OPTIONAL ] fill color (default:  light blue)
%
% See also SHARPLOEWNERMAJORANTP2, BETADOMAINP2.
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

%% setup
tol     = 1e-9;
m       = length(p);
Q       = orth( eye(m) - p*p');
mu_p    = p'*diag(Lambda)*p;
mu_hat  = max(eig(Q'*diag(Lambda)*Q));

b2_max  = 2*max(beta2vec)-min(beta2vec);
b1_max  = max(beta1vec);

if( ~exist('c', 'var') )
    c = [0,190,255]/255;
end

%% case 1: p is parallel to some eigenvector of A
if( idx_case == 1 )
    % replace beta1vec, beta2vec by a compact vector
    beta1vec = [ mu_hat, mu_hat, b1_max ];
    beta2vec = [ b2_max, mu_p-mu_hat, mu_p - b1_max ];
end

%% general plotting
h = figure;
hold on;
Fx = [ beta1vec, b1_max, beta1vec(1) ];
Fy = [ beta2vec,b2_max, b2_max];
fill( Fx, Fy, c, 'Edgecolor', 'none');

J           = find( abs( Lambda - Lambda(1) ) < tol * norm(Lambda) );
orthogonal  = ( norm(p(J)) < tol );

% if( abs( beta1vec(1) - mu_hat) < tol )
%     if( abs(mu_hat-Lambda(1)) < tol*norm(Lambda) && ~orthogonal && idx_case == 2 )
%         beta2vec(1) = b2_max;
%     else
%         beta2vec(2) = beta2vec(1);
%         beta2vec(1) = b2_max;
%         beta1vec(2) = beta1vec(1);
%     end
% end

if( idx_case == 1 )
    %% case 1: p is parallel to some eigenvector of A
    plot(beta1vec, beta2vec, '-', 'linewidth', 3, 'color', 'blue');
    xlabel('\beta_1');
    ylabel('\beta_2');
    title('case 1: p is contained in the eigenspace of A');
    a=axis();
else
    % if p orthogonal to the eigenspace of Lambda(1) ???
    
    %% case 2: p is inclined w.r.t. the eigensystem of A
    plot(beta1vec, beta2vec, '-', 'linewidth', 3, 'color', 'blue');
    
    hold on;
    %% special consideration:
    % p is inclined to the eigensystem of A AND p is orthogonal to the
    % Q-eigenvectors of mu_hat --> add point at (mu_hat, 'inf')
    if( orthogonal )
        % Case 
        plot([beta1vec(1), beta1vec], [b2_max, beta2vec], '-', 'linewidth', 3, 'color', 'blue');
        title('case 2: p is inclined to the eigensys. and orth. to Q_{max}');
    else
        title('case 2: p is inclined to the eigensystem of A');
    end
    
    xlabel('\beta_1');
    ylabel('\beta_2');
    a=axis();

end

plot([mu_hat, mu_hat, b1_max], [b2_max, mu_p-mu_hat, mu_p-b1_max],'--', 'color', [0.5,0.5,0.5],  'linewidth', 2 );


end
