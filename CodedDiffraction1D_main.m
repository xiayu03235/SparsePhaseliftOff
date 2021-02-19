%% Implementation of the Sparse Phaseliftoff algorithm proposed in the paper under coded diffraction patterns
%  ''Sparse phase retrieval via Phaseliftoff'' 
% by Y. Xia and Z. Q. Xu.

clc;
clear all;
close all;

%% Parameter setting
Params.n           = 50;                                    % signal dimension
Params.k           = 5;                                     % sparsity of signal
Params.L           = 3;                                     % number of diffraction patterns
Params.m           = Params.n * Params.L;                   % number of measurements
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
x                           = randn(Params.n, 1); 
loc                         = randperm(Params.n);
x(loc(Params.k + 1: end))   = 0;
x                           = x / norm(x);   % normalize the input

% measurement generation
A=[];
for L_num=1 : Params.L
    Masks = randsrc(Params.n, 1, [1i -1i 1 -1]);
 
    temp = rand(size(Masks));
    Masks = Masks .* ((temp <= 0.2) * sqrt(3) + (temp > 0.2) / sqrt(2));
    A=[A ; fft(eye(Params.n)) * diag(Masks)];
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