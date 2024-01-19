function [U, S, V] = tsvd(X, trank)


[n1, n2, n3] = size(X);
X = fft(X, [], 3);

min12 = min(n1, n2);
U = zeros(n1, min12, n3);
S = zeros(min12, min12, n3);
V = zeros(n2, min12, n3);

% i=1
[U(:, :, 1), S(:, :, 1), V(:, :, 1)] = svd(X(:, :, 1), 'econ');
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2:halfn3
    [U(:, :, i), S(:, :, i), V(:, :, i)] = svd(X(:, :, i), 'econ');
    U(:, :, n3+2-i) = conj(U(:, :, i));
    V(:, :, n3+2-i) = conj(V(:, :, i));
    S(:, :, n3+2-i) = S(:, :, i);
end
% if n3 is even
if mod(n3, 2) == 0
    i = halfn3 + 1;
    [U(:, :, i), S(:, :, i), V(:, :, i)] = svd(X(:, :, i), 'econ');
end

U = U(:, 1:trank, :);
V = V(:, 1:trank, :);
S = S(1:trank, 1:trank, :);

U = ifft(U, [], 3);
S = ifft(S, [], 3);
V = ifft(V, [], 3);
