function [A, L, S] = acc_GNBRL_solver(X, opts)

%% Parameters and defaults
if isfield(opts, 'maxit')     maxit = opts.maxit;   else maxit = 50;            end
if isfield(opts, 'tol')       tol = opts.tol;       else tol = 1e-2;            end
if isfield(opts, 'alpha')     alpha = opts.alpha;   else alpha = [1, 1, 1]/3;   end
if isfield(opts, 'Beta')      Beta = opts.Beta;     else Beta = 0.05;           end
if isfield(opts, 'rho')       rho = opts.rho;       else rho = 1.3;             end
if isfield(opts, 'lambda')    lambda = opts.lambda; else lambda = [0.5, 0.05];  end

psi = opts.psi; A = opts.A;

%% Data preprocessing and initialization
Nway = size(X); n1 = Nway(1); n2 = Nway(2); n3 = Nway(3); r = n2;
L = ones(r, n2, n3); hat_L = L; L0 = L;
hat_A = A; A0 = A;
S = X - A;

for u = 1:3 T{u} = zeros(n1, r, n3);
    beta(u) = 1e-4;
end
for u = 1:3
    C{u} = diff(A, 1, u);
end
C{1}(n1, :, :) = A(1, :, :) - A(end, :, :);
C{2}(:, r, :) = A(:, 1, :) - A(:, end, :);
C{3}(:, :, n3) = A(:, :, 1) - A(:, :, end);
epsilon = eps; l0 = 1; rw = 0.9999; gamma = 1.0001;
L_A0 = 1; L_L0 = 1;

%%
Eny_x = (abs(psf2otf([+1; -1], [n1, r, n3]))).^2;
Eny_y = (abs(psf2otf([+1, -1], [n1, r, n3]))).^2;
Eny_z = (abs(psf2otf([+1, -1], [r, n3, n1]))).^2;
Eny_z = permute(Eny_z, [3, 1, 2]);

fprintf('Iteration:     ');
for t = 1:maxit
    fprintf('\b\b\b\b\b%5i', t);

    S0 = S;

    %% update A
    L_A = tsn(L)^2 + epsilon;
    nabla_fA = tprod(tprod(hat_A, L)+S-X, tran(L));
    LAf = L_A * hat_A - nabla_fA;
    for u = 1:3 CT{u} = beta(u) * C{u} - T{u}; end
    diffHCT = diffH3(CT, [n1, r, n3]);
    A = real(ifftn(fftn(diffHCT+Beta*LAf)./(Beta * L_A + beta(1) * Eny_x + beta(2) * Eny_y + beta(3) * Eny_z)));

    %% 计算差分 nabla_A
    for u = 1:3
        nabla_A{u} = diff(A, 1, u);
    end
    nabla_A{1}(n1, :, :) = A(1, :, :) - A(end, :, :);
    nabla_A{2}(:, r, :) = A(:, 1, :) - A(:, end, :);
    nabla_A{3}(:, :, n3) = A(:, :, 1) - A(:, :, end);

    %% update C_u
    for u = 1:3
        nabla_ATu = nabla_A{u} + T{u} / beta(u);
        [C{u}, ~, ~] = prox_tnn_psi(nabla_ATu, alpha(u)/beta(u), psi);
    end

    %% update L
    nabla_fL = tprod(tran(A), tprod(A, hat_L)+S-X);
    L_L = gamma * tsn(A)^2 + epsilon;
    Lnabla_f = hat_L - nabla_fL / L_L;
    [L, ~, ~] = prox_tnn_psi(Lnabla_f, lambda(1)/L_L/Beta, psi);

    %% update S
    Z = X - tprod(A, L);
    for i = 1:n1
        for j = 1:n2
            z = Z(i, j, :);
            norm_z = norm(z(:));
            ss = z / norm_z * prox_psi(norm_z, lambda(2)/Beta, psi);
            S(i, j, :) = ss;
        end
    end

    %% update T_u and beta(u)
    for u = 1:3
        T{u} = T{u} + beta(u) * (nabla_A{u} - C{u});
        beta(u) = rho * beta(u);
    end
    Beta = min(1.5*Beta, 1e8);

    %% do extrapolation
    l = (1 + sqrt(1+4*l0^2)) / 2;
    w = (l0 - 1) / l; % extrapolation weight
%   w = 0;  %  no acceleration
    wA = min([w, rw * sqrt(L_A0/L_A)]);
    cgamma = max((gamma - 1)/2/gamma, 0.99^t);
    wL = min([w, cgamma * rw * sqrt(L_L0/L_L)]);
    hat_A = A + wA * (A - A0);
    hat_L = L + wL * (L - L0);

    % store old update
    l0 = l; A0 = A; L_A0 = L_A; L0 = L; L_L0 = L_L;

    %% judge whether converges
    vartheta2 = 0;
    for u = 1:3
        vartheta2 = max(vartheta2, norm(nabla_A{u}(:)-C{u}(:))/norm(C{u}(:)));
    end
    if vartheta2 < tol

        % Compute vartheta_A
        DT = diffH3(T, [n1, r, n3]); nabla_fA = tprod(tprod(A, L)+S-X, tran(L));
        vartheta_A = norm(DT(:)+Beta*nabla_fA(:)) / (1 + norm(DT(:))^2 + norm(Beta*nabla_fA(:))^2);

        % Compute vartheta_L
        nabla_fL = tprod(tran(A), tprod(A, L)+S-X);
        [eL, ~, ~] = prox_tnn_psi(L-Beta*nabla_fL, lambda(1), psi);
        vartheta_L = norm(L(:)-eL(:)) / (1 + norm(L(:))^2 + norm(Beta*nabla_fL(:))^2);

        % Compute vartheta_S
        nabla_fS = tprod(A, L) + S - X; eZ = S - Beta * nabla_fS;
        for i = 1:n1
            for j = 1:n2
                ez = eZ(i, j, :);
                norm_ez = norm(ez(:));
                ss = ez / norm_ez * prox_psi(norm_ez, lambda(2), psi);
                eS(i, j, :) = ss;
            end
        end
        vartheta_S = norm(S(:)-eS(:)) / (1 + norm(S(:))^2 + norm(Beta*nabla_fS(:))^2);

        % Compute vartheta_Cu
        vartheta_C = 0;
        for u = 1:3
            eC = prox_tnn_psi(C{u}+T{u}, alpha(u), psi);
            vartheta_C = max(vartheta_C, norm(C{u}(:)-eC(:))/(1 + norm(C{u}(:))^2 + norm(T{u}(:)))^2);
        end

        vartheta1 = [vartheta_A, vartheta_L, vartheta_S, vartheta_C];
        if max(vartheta1(:)) < 1e-1
            break;
        end
    end

    %    res = norm(S(:)-S0(:)) / norm(S0(:));

%         resa(t) = AUC(S, mask)

end
end