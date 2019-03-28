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
% Aug 16, 2018      initial creation
%

n_points    = 1000;     % number of points for the discrete value for B_2
n_interp    = 1000;     % number of interpolation points for the discrete value for B_2
n           = 10;       % dimension of the matrix
m           = 5;        % number of matrices to construct
tol         = 1e-8;
b2cutoff    = 1e4;     % very high number since otherwise,
%                         beta1-->mu_hat leads to numerical
%                         issues in the validation
b1_max      = n;

%% define one random vector p0
p0  = randn(n,1);
p0  = p0 / norm(p0);

%% generate a set of random matrices
A       = cell(m,1);
U       = cell(m,1);
l       = cell(m,1);
p       = cell(m,1);
a       = cell(m,1);
Lambda  = cell(m,1);
mu_hat  = cell(m,1);
mu_p    = cell(m,1);
idx_case= cell(m,1);

B2crit = []; % critical B2 domain for the set of all matrices
for i=1:m
    A{i}= randn(n+1,n+1);
    A{i}= A{i} + A{i}' + 0.05*n*eye(n+1); % make the matrix possibly indefinite
%     A{i}= A{i} + A{i}' - 0.25*n*eye(n+1); % make the matrix possibly indefinite
    %% split the matrix
    A2      =   A{i}( 1:n, 1:n );
    a{i}    =   A{i}( 1:n, n+1 );
    a0{i}   =   A{i}( n+1, n+1 );

    %% majorant for the upper n x n - block matrix
    [ U{i}, Lambda{i}, p{i}, mu_p{i}, mu_hat{i}, idx_case{i} ] = SharpLoewnerMajorantP2( A2, p0, tol );

    %% get B2(A2) --> update critical domain
    [b1, b2 ] = BetaDomainP2(  mu_hat{i}, b1_max, Lambda{i}, p{i}, idx_case{i}, n_points, b2cutoff );
    if( i == 1 )
        BETA    = [ b1; b2 ];
    else
        [b1, b2]= IntersectMajorants( b1, b2, BETA(1,:), BETA(2,:), n_interp );
        BETA    = [b1;b2];
    end
        
    %% transform a into l
    l{i}  = U{i}'*a{i};
end

%% optional - store data
% save('../data/paper_b4_set_optimization.mat');
%% optional - plot the admissible beta1/2 domain B2
% FillPlotB2(BETA);
% figure;
%% use data from the article:
load('../data/paper_b4_set_optimization.mat');

a0_max = max( [ a0{:} ] );



%% Semidefinite programming problem
% objective function: Frobenius norm squared
phi = @(beta) sqrt(n*beta(1)^2 + beta(2)^2 + beta(3)^2 + 2*beta(4)^2);

% compute discrete values of phi along the double sharp majorants
% computed earlier
%% initial guess:
% minimum over beta_1/beta_2 + delta_1/2 for beta_4=0, beta_3=beta_3^c
delta1  = 0.05;
delta2  = 0.05;
F       = zeros(size(BETA,2),1);
emin    = zeros(size(BETA,2),1);
DATA    = zeros(4,size(BETA,2));
DATA(1:2,:) = BETA(:,:);
for i=1:size(BETA,2)
    % initial guess for beta1/2: slightly inside of the critical two-parametric domain
    b1  = BETA(1,i) + delta1;
    b2  = BETA(2,i) + delta2;
    % then beta3 can be chosen based on the assertion beta4=0
    b3  = a0{1};
    for imat=1:m
        b3 = max( b3, a0{imat} + l{imat}' * (( diag(b1-Lambda{imat})+b2*p{imat}*p{imat}' )\l{imat} ));
    end
    b4  = 0;
    F(i) = phi([b1,b2,b3,b4]);
%     emin(i) = MinSetEigenvalue( [b1;b2;b3;b4], Lambda, p, l, a0 );
    DATA(1,i)=b1;
    DATA(2,i)=b2;
    DATA(3,i)=b3;
    DATA(4,i)=b4;
end
subplot(1,2,1);
plot(emin,'-', 'color', 'magenta', 'linewidth', 2);
% plot phi over beta1
subplot(1,2,2);
plot(DATA(1,:), F, '-', 'linewidth', 2, 'color', 'black');
% get optimum
[Fmin, I]   = min(F);
bmin        = DATA(:,I);
bmin0       = bmin + [ 0.1; 0.1; 0.1; 0];
Fmin0       = Fmin;
axis([DATA(1,1)-0.5 b1_max Fmin 8*Fmin]);
%% use constraint optimization to find (possibly sharp) majorants minimizing
% the objective function using numerical optimization

%% define **some** options for FMINCON
% Set the Display option to 'iter', the StepTolerance
% and the constraint tolerance (!!! important !!!)
options = optimoptions('fmincon');
options.Display                 = 'iter';
options.StepTolerance           = 1e-14;
options.ConstraintTolerance     = 1e-12;
options.MaxFunctionEvaluations  = 4000;
options.Algorithm               = 'active-set';
option.OptimalityTolerance      = 1e-12;
option.FunctionTolerance        = 1e-8;

%% define constraints (here: pass needed parameters to the function)
C4      = @(x,imat) [ diag(x(1)-Lambda{imat})+x(2)*p{imat}*p{imat}', x(4)*p{imat}-l{imat}; x(4)*p{imat}'-l{imat}', x(3)-a0{imat} ];
cond    = @(x) conditions( x,Lambda, p, l, a0);
phi     = @(x) Phi_Frobenius(x, n);
[b1_min, imin]  = max([mu_hat{:}]);
if( idx_case{imin} ~= 1 )
    b1_min = b1_min + 1e-6; % do not allow beta1 to be 'on' the lower limit
    % since we are in case 2 of the B2-Domain (p inclined)
end

ok      = 0;    
Fmin    = inf;
x0      = bmin0;
it      = 0;
while( ok == 0 )
    it = it + 1;
    [ bmin_loc, Fmin_loc, exitflag ] = fmincon( phi, x0, [], [], [], [], [b1_min;-inf;-inf;-inf], [], cond, options );
    if( exitflag > 0 )
        bmin = bmin_loc;
        Fmin = Fmin_loc;
        ok = 1;
    else
        if( max( cond( bmin_loc) ) > 0 )
            bmin_loc = x0;
        end
        
        test = 0;
        while( test == 0 )
            x0 = bmin_loc + [0; 0.2 * rand(2,1); 0];
            c = cond(x0);
            if( c(1) <= 0 && c(2) <= 0 )
                test = 1;
            end
        end
        fprintf('\n\n\n%% RESTARTING OPTIMIZATION\n');
    end
end
fprintf('\n\n**********************************************************************************************\n');
fprintf('* optimal value for double sharp Loewner Majorants   %16.8e\n', Fmin0);
fprintf('*          beta = [ %12.6e; %12.6e; %12.6e; %12.6e ];\n', ...
    bmin0(1), bmin0(2), bmin0(3), bmin0(4) );
emin = MinSetEigenvalue( bmin0, Lambda, p, l, a0);
fprintf('* smallest eigenvalue of the difference matrix: %12.5e\n', emin);
fprintf('**********************************************************************************************\n\n\n');
fprintf('**********************************************************************************************\n');
fprintf('* optimal value for all Loewner Majorants            %16.8e\n', Fmin);
fprintf('*          beta = [ %12.6e; %12.6e; %12.6e; %12.6e ];\n', ...
    bmin(1), bmin(2), bmin(3), bmin(4) );

emin = MinSetEigenvalue( bmin, Lambda, p, l, a0);
fprintf('* smallest eigenvalue of the difference matrix: %12.5e\n', emin);
fprintf('**********************************************************************************************\n\n\n');
fprintf('**********************************************************************************************\n');
fprintf('* change of minimum over double sharp majorants         %16.8e (relative: %10.5f)\n', Fmin-Fmin0, Fmin/Fmin0-1);
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
m   = length(Lambda);
c   = [-inf, -inf];
b2c = zeros(1, m);
for imat=1:m
    n           = p{imat}./(x(1)-Lambda{imat});
    b2c(imat)   =  - 1 / (n'*p{imat});

    z   = x(4)*p{imat} - l{imat};
    %% point on the boundary of the critical B2 domain!
    if( abs(x(2) - b2c(imat)) < tol )
        % assert x(2) = b2c
%         b4crit  = -b2c(imat) * (n'*l{imat});
        d       = abs(n'*z); % must be zero
        if( d > tol ) 
            %% beta4 is inadmissible --> penalty function
            fac  = 1e8;
            %% Attention: fac and 1/tol should not differ by more than
            % sqrt(accuracy_fmincon)
            c(2) = max( c(2), d * fac);
        else
            %% beta4 is admissible --> use pseudo-inverse
            c(2)= max( c(2), a0{imat} + sum( z.^2 ./ (x(1)-Lambda{imat}) ) - x(3) );
        end
    else
        if( x(2) < b2c(imat) )
            c(1)= max( c(1), b2c(imat)-x(2) );
            c(2)= max( c(2), 1 );
        else
            % db3 = z. C_2^-1 * z    (using Sherwood-Morrison-Woodbury)
            db3 = sum( z.^2./(x(1)-Lambda{imat}) ) - x(2)/(1 - x(2)/b2c(imat)) * (z'*n)^2;
            % (verified!)
            c(2)= max( c(2), a0{imat} + db3 - x(3) );
        end
    end
end

c(1) = max(  b2c(:) - x(2) );

%     z   = x(4)*p{imat} - l{imat};
% 
%     %% check if beta2=beta2_crit(beta1) --> C has a nullspace!
%     if( abs( x(2) + 1/(n'*p{imat}) ) < tol )
%         %% check if beta4 is admissible (i.e. such that a Loewner Majorant is
%         % still obtained)
%         
%         
%         if( d > tol )
%             %% beta4 is inadmissible --> penalty function
%             fac  = 1e12;
%             %% Attention: fac and 1/tol should not differ by more than
%             % sqrt(accuracy_fmincon)
%             c(2) = max( c(2), d * fac);
%             return;
%         else
%             %% check if beta3 is ok
%             zi  = pinv( diag(x(1)-Lambda{imat}) + x(2)*p{imat}*p{imat}') * z;
%             c(2)= max( c(2), a0{imat} + z'*zi - x(3) );
%             return;
%         end
%     else
%         if( x(2) < - 1/(n'*p{imat}) - tol )
%             c(1) = - 1/(n'*p{imat}) + tol - x(2);
%         else
%             %% check if constraint for beta3 is fulfilled
%             % (note: zi computes exactly since the matrix is positive definite)
%             z   = x(4)*p{imat} - l{imat};
%             zi  = ( diag(x(1)-Lambda{imat}) + x(2)*p{imat}*p{imat}') \ z;
%             c(2)= max( c(2), a0{imat} + z'*zi - x(3) );
%         end
%     end
% end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function emin = MinSetEigenvalue( x, Lambda, p, l, a0 )
m   = length(Lambda);
emin= zeros(m,1);
b1  = x(1);
b2  = x(2);
b3  = x(3);
b4  = x(4);
for imat=1:m
    emin(imat) = min(eig( [ diag(b1-Lambda{imat})+b2*p{imat}*p{imat}', b4*p{imat}-l{imat}; b4*p{imat}'-l{imat}', b3-a0{imat} ] ) );
end
% emin
emin=min(emin);
end
