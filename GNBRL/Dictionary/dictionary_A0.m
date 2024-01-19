function L = dictionary_A0(X, psi)
% Solve the Tensor Robust Principal Component Analysis based on Tensor Nuclear Norm problem by ADMM
%
% min_{L,S} ||L||_*+lambda*||S||_1, s.t. X=L+S
%
% ---------------------------------------------
% Input:
%       X       -    d1*d2*d3 tensor
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       L       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual
%       iter    -    number of iterations
%
% version 1.0 - 19/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
%

tol = 1e-8;
max_iter = 100;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
lambda = 0.005;  %0.01,0.01,0.01,0.01,0.005,0.005
DEBUG = 0;
[n1,n2,n3] = size(X);
if ~exist('opts', 'var')          opts = [];                 end
if isfield(opts, 'tol');          tol = opts.tol;            end
if isfield(opts, 'max_iter');     max_iter = opts.max_iter;  end
if isfield(opts, 'rho');          rho = opts.rho;            end
if isfield(opts, 'mu');           mu = opts.mu;              end
if isfield(opts, 'max_mu');       max_mu = opts.max_mu;      end
if isfield(opts, 'lambda');       lambda = opts.lambda;      end
if isfield(opts, 'DEBUG');        DEBUG = opts.DEBUG;        end

dim = size(X);
L = zeros(dim);
S = L;
Y = L;
L = X;

for iter = 1:max_iter
    Lk = L;
    Sk = S;
    
    % update L
%    [L, tnnL, trank] = prox_tnn_psi(-S+X-Y/mu, 1/mu, psi);
%    [L, tnnL, trank] = prox_tnn_wg(-S+X-Y/mu, 1/mu);
     [L, tnnL, trank] = prox_tnn_A0(-S+X-Y/mu, 1/mu);
    
    % update S
    Z = -L+X-Y/mu;
    for i = 1:n1
        for j = 1:n2
            z = Z(i, j, :);
            norm_z = norm(z(:));
            ss = z / norm_z * prox_CapLog(norm_z, lambda/mu, 0.1, 0.1);
            S(i, j, :) = ss;
        end
    end

    dY = L + S - X;
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([chgL, chgS, max(abs(dY(:)))]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            %             obj = tnnL+lambda*norm(S(:),1);
            %             err = norm(dY(:));
            %             disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
            %                     ', obj=' num2str(obj) ', err=' num2str(err)]);
            sparsityhat = length(find(S ~= 0));
            fprintf('iter = %d, obj = %.3f, err = %f, mu=%.2f, rankL = %d, sparsity = %d\n' ...
                , iter, tnnL+lambda*norm(S(:), 1), chg, mu, trank, sparsityhat);
        end
    end

    if chg < tol
        break;
    end
    Y = Y + mu * dY;
    mu = min(rho*mu, max_mu);
end

end
