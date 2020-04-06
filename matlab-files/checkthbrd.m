Q = dlmread('outQ', '\t', 1,0);
D = dlmread('outD', '\t', 1,0);
E = dlmread('outE', '\t', 1,0);

T = diag(D) + diag(E,-1) + diag(E,1);
A = diag(-2*ones(25,1)) + diag(ones(24,1),1) + diag(ones(24,1),-1) + diag(ones(20,1),5) + diag(ones(20,1),-5);
%A = diag(-2*ones(4,1)) + diag(ones(3,1),1) + diag(ones(3,1),-1) + diag(ones(2,1),2) + diag(ones(2,1),-2);

fprintf('ZHBTRD error: %f\n', norm(Q' * A * Q - T));

EV = dlmread('outEV', '\t', 1,0);
EE = dlmread('outEE', '\t', 1,0);

[testV, testD] = eig(A);
fprintf('ZSTEQR eigenvalue error: %f\n', norm(testD - diag(EE)));
fprintf('ZSTEQR eigenvector error: %f\n', norm(abs(testV) - abs(EV)));
