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
%



%% set up the matrix A
A       = diag( [ 2.5, 1, 0.5 ] );
tol     = 1e-8;
b2cutoff= 100;
n_points= 400;
b1_max  = 30;
%% TEST CASE 1:
fprintf('%% ---------------------------------------------------\n');
fprintf('%% Test Case 2: p = e_1\n');
p0      = [ 0.6; 0.8; 0. ];

[ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A, p0, tol );
[ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, b1_max, Lambda, p, idx_case, n_points, b2cutoff );
MajorantInfoP2( Lambda, p, idx_case );
PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case );

%% selected results:
BETA = [    1.98,   25.48;
            2.0,    12.5;
            2.1,    3.1429;
            2.25,   1.0776;
            2.5,    0;
            3,     -0.9615;
            4,     -2.2059;
            8,     -6.3742; ];
hold on;
plot( BETA(:,1), BETA(:,2), 'o', 'color', 'black', 'markerfacecolor', 'black');

%% TEST CASE 2:
fprintf('\n%% ---------------------------------------------------\n');
fprintf('%% Test Case 2: p = e_1\n');
p0  = [ 1.0; 0.; 0. ];

[ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A, p0, tol );
[ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, b1_max, Lambda, p, idx_case, n_points, b2cutoff );
MajorantInfoP2( Lambda, p, idx_case );
PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case );

%% selected results:
BETA = [    1.0,  1.5;
            1.5,  1.0;
            2.5,  0.0;
            4.0, -1.5;
            8.0, -5.5; ];
plot( BETA(:,1), BETA(:,2), 'o', 'color', 'black', 'markerfacecolor', 'black');

%% TEST CASE 3:
fprintf('\n%% ---------------------------------------------------\n');
fprintf('%% Test Case 3: p = e_2 \n');
p0  = [ 0.0; 1.0; 0. ];

[ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A, p0, tol );
[ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, b1_max, Lambda, p, idx_case, n_points, b2cutoff );
MajorantInfoP2( Lambda, p, idx_case );
PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case );


%% selected results:
BETA = [    2.5, -1.5;
            3.0, -2.0;
            4.0, -3.0;
            8.0, -7.0; ];
plot( BETA(:,1), BETA(:,2), 'o', 'color', 'black', 'markerfacecolor', 'black');

%% TEST CASE 4:
fprintf('\n%% ---------------------------------------------------\n');
fprintf('%% Test Case 4: p = sqrt(0.5) * ( e_2 + e_3 )\n');
p0      = [ 0.; 1; 1 ]/sqrt(2);

[ U, Lambda, p, mu_p, mu_hat, idx_case ] = SharpLoewnerMajorantP2( A, p0, tol );
[ beta1vec, beta2vec ] = BetaDomainP2( mu_hat, 8, Lambda, p, idx_case, n_points, b2cutoff );
MajorantInfoP2( Lambda, p, idx_case );
PlotCriticalB2( beta1vec, beta2vec, Lambda, p, idx_case );
