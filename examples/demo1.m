clc;
clear all;
close all;
addpath(genpath('..'))

%% Examples for Lemma 1,2,3

disp('*********************************************************')
disp('Examples for Lemma 1,2,3 (figures 1,2,3)')

% Matrix A
A = diag([2.5,1,0.5]);

% Vector p0 = p1
p0 = [1,0,0]';
figure;
B2Plot(A,p0,1,1,[0,5],[-6,10]);
legend('Boundary of B2 for p1','Necessary condition','Necessary condition')

% Vector p0 = p2
p0 = [0,1,1]';
p0 = p0/norm(p0);
figure;
subplot(1,2,1);
B2Plot(A,p0,1,1,[0,5],[-6,10]);
legend('Boundary of B2 for p2','Necessary condition','Necessary condition')
subplot(1,2,2);
B2Plot(A,p0,0,1,[2.45,2.55],[-1.9,-1.5]);
legend('Boundary of B2 for p2','Necessary condition','Necessary condition')

% Vector p0 = p3
p0 = [3,4,0]';
p0 = p0/norm(p0);
figure;
B2Plot(A,p0,1,1,[0,5],[-6,10]);
legend('Boundary of B2 for p3','Necessary condition','Necessary condition')

%% Further examples for Lemma 2

disp('*********************************************************')
disp('Further examples for Lemma 2 (figure 4)')

figure;

% muQ = L1 and muQ <= b1
A = diag([3,3,2,5/10]);
p0 = [1,0,0,1]';
p0 = p0/norm(p0);
subplot(1,3,1);
B2Plot(A,p0,1,1,[0,5],[-6,10]);
legend('Boundary of B2','Necessary condition','Necessary condition')
title('Lemma 2 - mu_Q = L1 and mu_Q <= b1')

% muQ = L2 and muQ <= b1
A = diag([3,2,2,5/10]);
p0 = [1,0,0,1]';
p0 = p0/norm(p0);
subplot(1,3,2);
B2Plot(A,p0,1,1,[0,5],[-6,10]);
legend('Boundary of B2','Necessary condition','Necessary condition')
title('Lemma 2 - mu_Q = L2 and mu_Q <= b1')

% muQ = L2 and muQ<b1
A = diag([3,2,2,1]);
p0 = [1,0,0,1]';
p0 = p0/norm(p0);
subplot(1,3,3);
B2Plot(A,p0,1,1,[0,5],[-6,10]);
legend('Boundary of B2','Necessary condition','Necessary condition')
title('Lemma 2 - mu_Q = L2 and mu_Q < b1')