function [X_out] = ADMM_sub(A, y, Y, lambda, mu, n, m) 
%%Input Parameters:
% A                 the transpose of measurement matrix
% y                 the obeservation
% Y                 the inner variable in SparsePhaseliftoff
% lambda, mu        the parameters in the model
% m                 the dimension of measurements
% n                 the dimension of signal 

%%Output:
% X_out             the output matrix


%% ADMM initialization
X1       = zeros(n,n);
X2       = zeros(n,n);
X3       = zeros(n,n);

Y1       = zeros(n,n);
Y2       = zeros(n,n);

delta    = 1;        % parameter in ADMM
iter_max = 3000;     % max number of iterations;

%% Loop iteraion
for iter = 1:iter_max
    X1_temp             = A * diag(y) * A' - Y1 - lambda * (eye(n) - Y) + delta * X3;
    X1                  = 1 / delta * (X1_temp - A * diag(inv(abs(A' * A) .^ 2 + delta * eye(m)) * diag(A' * X1_temp * A)) * A');
    X2                  = (X3 -1 / delta * Y2) .* max(1 - (mu / delta)./abs(X3 -1 / delta * Y2), 0);
    X3_temp             = 1 / 2 * (X1 + X2) + 1 / (2 * delta) * (Y1 + Y2);
    [U , Sigma, V]      = eig(X3_temp);
    X3                  = U * max(real(Sigma), 0) * U';
    Y1                  = Y1 + delta * (X1 - X3);
    Y2                  = Y2 + delta * (X2 - X3);  
end
X_out = X3;
end
