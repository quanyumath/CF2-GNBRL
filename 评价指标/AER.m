function f = AER(E, mask)

% AER 值越大越好 （The larger the AER value, the better the detection performance）

r_max = max(E(:)); num = size(mask, 1) * size(mask, 2);
taus = linspace(0, r_max, 5000);
mask_reshape = reshape(mask, 1, num);
anomaly_map = logical(double(mask_reshape) > 0);
normal_map = logical(double(mask_reshape) == 0);
for index2 = 1:length(taus)
    tau = taus(index2);
    anomaly_map_rx = (E(:) > tau)';
    PF0(index2) = sum(anomaly_map_rx & normal_map) / sum(normal_map);
    PD0(index2) = sum(anomaly_map_rx & anomaly_map) / sum(anomaly_map);
end
norm_taus = normalize(taus,'range');
%norm_taus = taus;
APFA = sum((norm_taus(2:end) - norm_taus(1:end-1)).*(PF0(2:end) + PF0(1:end-1))/2);
APD = sum((norm_taus(2:end) - norm_taus(1:end-1)).*(PD0(2:end) + PD0(1:end-1))/2);
f = (1 - APFA) / (1 - APD);
end








