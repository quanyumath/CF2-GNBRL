function [X, tnn, trank] = prox_ttnn_psi(U, Y, rho, psi)


[n1, n2, n3] = size(Y);
max12 = max(n1, n2);
X = zeros(n1, n2, n3);
O = tenmat(Y, [3]); %square norm
SS = O.data;
U3 = Unfold(U, [n1, n2, n3], 3);
[UU, ~, ~] = svd(U3, 'econ');
Y = UU' * SS;
Y = tensor(tenmat(Y, O.rdims, O.cdims, O.tsize));
Y = Y.data;
tnn = 0;
trank = 0;

for i = 1:n3
    [U, S, V] = svd(Y(:, :, i), 'econ');
    S = diag(S);
    S = prox_psi(S, rho, psi);
    S = max(S-rho, 0);
    r = length(find(S ~= 0));
    S = S(1:r);
    X(:, :, i) = U(:, 1:r) * diag(S) * V(:, 1:r)';
    tnn = tnn + sum(S) * 2;
    trank = max(trank, r);
end

O = tenmat(X, [3]); %square norm
SS = O.data;
Y = UU * SS;
Y = tensor(tenmat(Y, O.rdims, O.cdims, O.tsize));
X = Y.data;

tnn = tnn / n3;
