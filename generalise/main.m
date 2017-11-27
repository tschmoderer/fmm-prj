clearvars;
clc;

N  = 5;
dx = 1/N;

n = 1;
t = 0;

theta = ones(n,n);

% [X Y] = meshgrid(-5:0.1:5, -5:0.1:5);
% Z = X.^2+Y.^2-1;
% %Z = peaks(X,Y);
% Omega = contourf(X,Y,Z,[0,0]);
% 


T = Inf*ones(N+3,N+3);
Theta = -1*ones(N+3,N+3);

%pts acceptes 
ap = zeros(N+1,N+1);
ap(1:2,:) = 1;

AP = zeros(N+3);
AP(2:end-1,2:end-1) = ap;
AP = logical(AP)

T(AP) = 0
Theta(AP) = 1
t     = 0;
