%% Implementation of the Sparse Phaseliftoff algorithm proposed in the paper under Gaussian measurements
%  ''Sparse phase retrieval via Phaseliftoff'' 
% by Y. Xia and Z. Q. Xu.

clc;
clear all;
close all;

%% Parameter setting
Params.n           = 50;                                    % signal dimension
Params.k           = 5;                                     % sparsity of signal
Params.m           = 60;                                    % number of measurements
Params.cplx_flag   = 0;                                     % real: cplx_flag = 0;  complex: cplx_flag =1
Params.noise_flag  = 0;                                     % noiseless: noise_flag = 0; noisy: noise_flag = 1
Params.T           = 30;                                    % number of iterations
Params.mu          = 1e-3;                                  % parameter mu in the model
Params.lambda      = Params.mu * Params.k / (sqrt(2) - 1);  % parameter mu in the model
Params.snr         = 0;                                     % noise level in dB when Params.noise_flag  = 1
Params.iter_max    = 30;                                    % maximum number of iterations
Params.tol         = 1e-3;                                  % stop tolerance for the algorithm
display(Params);

%% Make signal and observations
% sparse signal generation
x                           = randn(Params.n, 1) + Params.cplx_flag * 1i * randn(Params.n, 1); 
loc                         = randperm(Params.n);
x(loc(Params.k + 1: end))   = 0;
x                           = x / norm(x);   % normalize the input

% measurement generation
if Params.cplx_flag == 0
    A  = randn(Params.m, Params.n); % real measurements
else
    A  = (randn(Params.m, Params.n) + 1i * randn(Params.m, Params.n)) / sqrt(2); % complex measurements
end

% observation generation
if Params.noise_flag == 0
    y             = abs(A * x) .^ 2; % noiseless measurements
else
    y0            = abs(A * x) .^ 2;
    y             = awgn(y0, Params.snr); % noisy measurements
    y(y<0)        = 0;
    Params.mu     = max(0.5 * norm(y - y0, 'fro'), 0.001);
    Params.lambda = Params.mu * Params.k / (sqrt(2) - 1);  % reset the parameters in the model
end

%% run Sparse Phaseliftoff algorithm
[Relerrs, z] = SparsePhaseliftOff(A, x, y, Params.lambda, Params.mu, Params.m, Params.n, Params.iter_max, Params.tol); 
disp('----------Sparse Phaseliftoff done!----------');

%% plot the relative error of Sparse Phaseliftoff
if length(Relerrs)>1
    figure,
    semilogy([1:length(Relerrs)], Relerrs)
    xlabel('Iteration'), ylabel('Relative Error'), ...
        title('Sparse Phaseliftoff: Relerr vs. itercount')
    grid
end