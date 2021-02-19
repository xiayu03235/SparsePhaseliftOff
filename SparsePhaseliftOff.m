function [Relerrs, z]=SparsePhaseliftOff(A, x, y, lambda, mu, m, n, iter_max,tol)
%%Input Parameters:
% A                 the measurement matrix
% x                 the input signal
% y                 the obeservation
% lambda, mu        the parameters in the model
% m                 the dimension of measurements
% n                 the dimension of signal
% iter_max          the maximum number of iterations
% tol               the stop tolerance 

%%Output:
% Relerrs        the error in each iterations
% z                 the output signal by SparsePhaseliftoff

%% Initialization

Relerrs    = [];
Y          = zeros(n,n);

%% Loop

for iter = 1:iter_max
    X             = ADMM_sub(A', y, Y, lambda, mu, n, m);
    Y             = X / norm(X, 'fro');
    [U, Sigma, V] = svd(X);        % since X is positive semidefinite, svd is the same as eigenvalue decomposition
    z             = U(:,1) * sqrt(Sigma(1,1));
    Relerr        = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro');
    Relerrs       = [Relerrs; Relerr];
    fprintf('Relative error in %d iteration is: %f\n', iter, Relerr);
    if Relerr < tol
        break;
    end
end  
end

