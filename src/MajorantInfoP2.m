function MajorantInfoP2( Lambda, p, idx_case )
%% MAJORANTINFOP2 displays some information on the Loewener majorant
% determined by SHARPLOEWNERMAJORANTP2.
%
% INPUTS     (remark: all inputs are outputs of SHARPLOEWNERMAJORANTP2)
% Lambda     vector containing the eigenvalues of A (sorted)
% p          vector p0 in the A-eigensystem
% idx_case   either 1 or two
%
%
%
% See also SHARPLOEWNERMAJORANTP2, BETADOMAINP2, PLOTCRITICALB2
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
% https://github.com/EMMA-Group/LoewnerMajorant
%
% CHANGELOG
% Aug 08, 2018      initial creation

%% setup
tol     = 1e-9;
m       = length(p);
Q       = orth( eye(m) - p*p');
mu_p    = p'*diag(Lambda)*p;
mu_hat  = max(eig(Q'*diag(Lambda)*Q));

%% display some info on the console
fprintf('%% ---------------------------------------------------\n');
fprintf('%% SHARP LOEWNER MAJORANT OF A SYMMETRIC MATRIX A\n');
fprintf('%% ---------------------------------------------------\n');
fprintf('%% dimension ......... %4i\n', m);
fprintf('%% mu_p = p''*A*p ..... %15.6e\n', mu_p);
fprintf('%% mu_hat ............ %15.6e\n', mu_hat);
fprintf('%% ---------------------------------------------------\n');
if( idx_case == 1 )
    fprintf('%% CASE 1: p is eigenvector of A\n');
fprintf('%% ---------------------------------------------------\n');
    fprintf('%% exact linear relation for beta_2\n');
    fprintf('%% beta_1    >=  mu_hat         =  %12.5e\n', mu_hat );
    fprintf('%% beta_2    >=  mu_p - beta_1  =  %12.5e - beta_1\n', mu_p );
    fprintf('%% beta_1    =   mu_hat  ==>  beta_2  >=  mu_hat  is sharp\n');
else
    fprintf('%% CASE 2: p is inclined to the eigensystem of A\n');
    fprintf('%% ---------------------------------------------------\n');
    fprintf('%% asymptotic linear relation for beta_2\n');
    fprintf('%% for beta_1 --> mu_hat  ==>  beta_2  -->  inf\n' );
    fprintf('%% for beta_1 --> inf     ==>  beta_2  -->  mu_p - beta_1\n' );
    fprintf('%% [ general nonlinear relation ]\n');
end
fprintf('%% ---------------------------------------------------\n');

end
