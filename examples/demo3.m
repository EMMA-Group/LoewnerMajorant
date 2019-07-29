clc;
clear all;
close all;
addpath(genpath('..'))

%% Matrix and vector data

disp('----------------------------------')
disp('Pre-computations')

% Matrix set and vector
% A = {diag([2,0.5,-3]),[[-4,14,-16];[14,5,2];[-16,2,-10]]/9};
A = [[-4,14,-16];[14,5,2];[-16,2,-10]]/9;
p0 = [1,0,0]';
n = length(p0);

% Discretization
[b1,b2] = meshgrid(0:0.1:10,-10:0.1:6);

% Local precomputed functions for determination of minimum
[idx,Lv,p,mu_P,mu_Q,c_Q] = B2Preparations(A,p0,1);
Dinv = @(b1) diag(1./(b1 - Lv));
Dinv2 = @(b1) diag(1./(b1 - Lv).^2);
b2min = @(b1) B2b2min(b1,idx,Lv,p,mu_P,mu_Q,c_Q);
db2min = @(b1) -(p'*Dinv2(b1)*p)/(p'*Dinv(b1)*p)^2;
tangent = @(b1) [1,db2min(b1)]';
eq = @(b1,dphi) dphi(b1,b2min(b1))'*tangent(b1);

%% Plots

fig = figure;
colormap((0.4:0.001:1)'*[1,1,1]);

%% Minimization for linear \varphi

disp('----------------------------------')
disp('Minimization for linear \varphi')

subplot(1,2,1);
daspect([1 1 1]);
hold on;

% Trace
phi = n*b1 + b2;

% Contour plot
contourf(b1,b2,phi,'LineStyle','none');

% Boundaries
B2Plot(A,p0,0,0,[0,10],[-10,6],200,'-','black');
title('$\varphi = 3 \beta_1 + \beta_2$','Interpreter','latex');

% Determine global minimum
dphi = @(b1,b2) [n,1]';
f = @(b1) eq(b1,dphi);
minb1 = fsolve(f,Lv(1)+0.001)
minb2 = b2min(minb1);
minphi = n*minb1 + minb2

% Plot minimum
scatter(minb1,minb2,'filled','black');
contour(b1,b2,phi,[minphi,minphi],'linewidth',2,'color','black','linestyle','--');

% Plot gradient and tangent
g = dphi(minb1,minb2);
t = tangent(minb1);
quiver(minb1,minb2,g(1),g(2),0,'color','black','linewidth',2,'maxheadsize',0.5);
quiver(minb1,minb2,t(1),t(2),0,'color','black','linewidth',2,'maxheadsize',0.5);

% Legend
colorbar;

hold off;

%% Minimization for nonlinear \varphi

disp('----------------------------------')
disp('Minimization for nonlinear \varphi')

subplot(1,2,2);
daspect([1 1 1]);
hold on;

% Trace
phi = sqrt(n*b1.^2 + b2.^2);

% Contour plot
contourf(b1,b2,phi,'LineStyle','none');

% Boundaries
B2Plot(A,p0,0,0,[0,10],[-10,6],200,'-','black');
title('$\varphi = \sqrt{3 \beta_1^2 + \beta_2^2}$','Interpreter','latex');

% Determine global minimum
dphi = @(b1,b2) 1/(2*sqrt(n*b1^2 + b2^2))*[n*2*b1,2*b2]';
f = @(b1) eq(b1,dphi);
minb1 = fsolve(f,Lv(1)+0.001)
minb2 = b2min(minb1);
minphi = sqrt(n*minb1^2 + minb2^2)

% Plot minimum
scatter(minb1,minb2,'filled','black');
contour(b1,b2,phi,[minphi,minphi],'linewidth',2,'color','black','linestyle','--');

% Plot gradient and tangent
g = dphi(minb1,minb2);
t = tangent(minb1);
quiver(minb1,minb2,g(1),g(2),0,'color','black','linewidth',2,'maxheadsize',0.9);
quiver(minb1,minb2,t(1),t(2),0,'color','black','linewidth',2,'maxheadsize',0.5);

% Legend
colorbar;

hold off;

%% Export

print('B2_examples','-dpng','-r600')