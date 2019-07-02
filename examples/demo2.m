clc;
clear all;
close all;
addpath(genpath('..'))

%% Example for matrix set with 2 matrices

disp('****************************************************');
disp('Example for matrix set with 2 matrices');

% Matrix set and vector
A = {diag([2,0.5,-3]),[[-4,14,-16];[14,5,2];[-16,2,-10]]/9};
N = length(A);
p0 = [1,0,0]';

% Plot individual B2 regions for b1 in [0,4]
figure;
subplot(1,2,1);
hold on;
for i=1:N
    B2Plot(A{i},p0,1,0,[0,4],[-6,10]);
end
legend('B2 for A1','B2 for A2');
title('Plot for each matrix')
hold off;

% Plot boundary of B2 for matrix set A for b1 in [0,4]
subplot(1,2,2);
B2Plot(A,p0,0,0,[0,4],[-6,10]);
title('Plot for matrix set')

%% Example for matrix set with N matrices in R^nxn

disp('****************************************************');
disp('Example for matrix set with N matrices in R^nxn');

% Set choose = 1 to load data and choose = 2 for random data generation
choose = 1;
if choose==1
    % Load prepared data
    load('data/A2set_p0.mat');
    N = length(A);
    n = length(p0);
    fprintf('Number of matrices in set: \n\t%i\nDimension of p0: \n\t%i',N,n);
end
if choose==2
    % Generate matrix set and vector p0
    N = 8;
    n = 100;
    fprintf('Number of matrices in set: \n\t%i\nDimension of p0: \n\t%i',N,n);
    A = cell(N,1);
    for i=1:N
        temp = rand(n);
        A{i} = (temp+temp')/2-n/35*rand()*eye(n);
    end
    p0 = rand(n,1);
    p0 = p0/norm(p0);
end

% Plot individual B2 regions
figure;
subplot(1,2,1);
hold on;
for i=1:N
    if choose==1
        B2Plot(A{i},p0,1,0,[0,20],[-20,30]);
    end
    if choose==2
        B2Plot(A{i},p0,1,0,[0,60],[-50,200]);
    end
end
title('Plot for each matrix')
hold off;

% Plot boundary of B2 for matrix set A
subplot(1,2,2);
if choose==1
    B2Plot(A,p0,0,0,[0,20],[-20,30]);
end
if choose==2
    B2Plot(A,p0,0,0,[0,60],[-50,200]);
end
title('Plot for matrix set')