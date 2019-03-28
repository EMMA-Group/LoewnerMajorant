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
% - example the majorization of a set of random symmetric matrices
% - illustrates how to access the different tools including
%   a discrete set majorization
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
%

%% random testing
d   = 100;          % size of the matrix
N   = 10;            % number of matrices to be considered

tol         = 1e-8;     % non-dimensional tolerance for colinearity detection
n_points    = 5000;     % number of points 
n_interp    =10000;     % number of points for the (fine grained) interpolation
b1_max      =  100;     % beta_1_max for the plotting and set majorization
b2_cutoff   =  300;
%% random direction p0
p0  = randn(d, 1);
p0  = p0/norm(p0);

%% BETA, BETA_c will contain the critical values
BETA = [];
figure;
hold on;
grid on;

%% define colors for the plotting
c   = cell(5,1);
c{1}= 'black';
c{2}= 'green';
c{3}= 'blue';
c{4}= 'yellow';
c{5}= 'magenta';

%% loop over N matrices and draw intersection
b1_min_global   = 0;
A               = cell(N,1);
Lambda          = cell(N,1);
p               = cell(N,1);
idx_case        = cell(N,1);
for i_mat=1:N
    fprintf('%% processing matrix  %5i out of %5i\n', i_mat, N );
    A{i_mat}    = randn(d,d);
    A{i_mat}    = A{i_mat} + A{i_mat}' + eye(d)*0.1;

    [ U, Lambda{i_mat}, p{i_mat}, mu_p, mu_hat, idx_case{i_mat} ] = SharpLoewnerMajorantP2( A{i_mat}, p0, tol );
    if( i_mat == 1 )
        b1_min_global = mu_hat;
    else
        b1_min_global = max( mu_hat, b1_min_global );
    end
    [ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, b1_max, Lambda{i_mat}, p{i_mat}, idx_case{i_mat}, n_points, b2_cutoff );
    plot(beta1vec, beta2vec, '-', 'color', c{mod(i_mat,5)+1}, 'linewidth', 2);
    if( i_mat == 1 )
        BETA = [ beta1vec; beta2vec ];
    else
        [b1_new, b2_new ] = IntersectMajorants( beta1vec, beta2vec, BETA(1,:), BETA(2,:), n_interp );
        BETA = zeros( 2, size(b1_new, 2) );
        BETA(1,:) = b1_new;
        BETA(2,:) = b2_new;
    end
end

%% re-evaluate the beta's on the new grid
%  (for file output and further processing)
allBETA= zeros( n_interp, N+2);
allBETA(:,1) = BETA(1,:);
allBETA(:,2) = BETA(2,:);
for i_mat=1:N
    [ ~, allBETA(:,i_mat+2) ] = BetaDomainP2( b1_min_global, b1_max, Lambda{i_mat}, p{i_mat}, idx_case{i_mat}, n_interp, b2_cutoff );
end
dlmwrite('all_beta_crit.dat', allBETA, 'precision', '%20.12e', 'delimiter', ' ');
A_all = zeros(N*d,d);
ct=1;
for i=1:N
    A_all(ct:(ct+d-1),:) = A{i};
    ct=ct+d;
end
dlmwrite('all_A.dat', A_all, 'precision','%20.12e', 'delimiter', ' ');
FillPlotB2(BETA);

% plot( BETA(1,:), BETA(2,:), '--', 'color', 'red', 'Linewidth', 4);
% grid on;
% axis( [ b1_min_global, b1_max, min(BETA(2,:)), max(BETA(2,:)) ] );
