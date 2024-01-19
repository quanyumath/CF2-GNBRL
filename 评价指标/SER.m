function f = SER(E, mask)

% SER 值越小越好 （The smaller the value of the SER, the better the
% anomaly detection performance of the algorithm）

[n1, n2] = size(mask);
Index_An = find(mask == 1); Index_Bg = find(mask == 0);
f = sum((E(Index_An) - 1).^2) + sum((E(Index_Bg) - 0).^2);
f = f / n1 / n2 * 100;

end