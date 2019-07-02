clc;
clear all;
close all;
addpath(genpath('..'))

%% Matrix and vector data

% Matrix set and vector
% A = {diag([2,0.5,-3]),[[-4,14,-16];[14,5,2];[-16,2,-10]]/9};
A = [[-4,14,-16];[14,5,2];[-16,2,-10]]/9;
p0 = [1,0,0]';
n = length(p0);

% Discretization
[b1,b2] = meshgrid(0:0.1:4,-6:0.1:10);

%% Contour plot of trace

figure;
subplot(1,2,1);
hold on;

% Trace
phi = n*b1 + b2;

% Contour plot
colormap((0:0.02:1)'*[1,1,1]);
contourf(b1,b2,phi);

% Boundaries
B2Plot(A,p0,0,0,[0,4],[-6,10]);
title('phi = 3 b1 + b2');

hold off;

%% Contour plot of Frobenius norm

subplot(1,2,2);
hold on;

% Trace
phi = sqrt(n*b1.^2 + b2.^2);

% % Contour plot
colormap((0:0.02:1)'*[1,1,1]);
contourf(b1,b2,phi);

% Boundaries
B2Plot(A,p0,0,0,[0,4],[-6,10]);
title('phi = sqrt(3 b1^2 + b2^2)');