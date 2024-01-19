function ROC(Methods_E, mask, xlab)

% The ROC describes the relationship between the x-axis FPR and
% y-axis TPR under different thresholds. The closer the curve
% is to the upper-left corner of the coordinate axis, the better
% the performance is.

% 说明：
% Methods_E---各种方法得到的异常矩阵
% DataTest----原始 HSIs 图像
% xlab--------各种方法的名称

Color{1} = [0, 0, 0]/255; Color{2} = [0, 255, 0]/255; Color{3} = [255, 0, 255]/255;
Color{4} = [128, 128, 255]/255; Color{5} = [255, 255, 0]/255; Color{6} = [0, 0, 255]/255;
Color{7} = [0, 255, 255]/255; Color{8} = [255, 204, 255]/255; Color{9} = [0, 64, 0]/255;
Color{10} = [255, 0, 0]/255;
% Color = cbrewer('seq', 'Greens', 10, 'linear');
for k = 1:length(Methods_E)
    E = Methods_E{k};
    if size(E, 3) > 1
        E = sqrt(sum(E.^2, 3));
        E = (E - min(E(:))) / (max(E(:)) - min(E(:)));
    end
    [PF0, PD0] = PFPD(E, mask);
    semilogx(PF0, PD0, 'Color', Color{k}, 'LineWidth', 3); hold on; % 或者将 semilogx 改成 plot
end
xlabel('False alarm rate', 'fontsize', 19);
ylabel('Probability of detection', 'fontsize', 19);
%legend(xlab);
axis([0, 1, 0, 1]);
hold off;

end

function [PF, PD] = PFPD(E, mask)
if size(E, 3) > 1
    E = sqrt(sum(E.^2, 3));
    E = (E - min(E(:))) / (max(E(:)) - min(E(:)));
end
r_max = max(E(:)); num = size(mask, 1) * size(mask, 2);
taus = linspace(0, r_max, 5000);
mask_reshape = reshape(mask, 1, num);
anomaly_map = logical(double(mask_reshape) > 0);
normal_map = logical(double(mask_reshape) == 0);
for index2 = 1:length(taus)
    tau = taus(index2);
    anomaly_map_rx = (E(:) > tau)';
    PF(index2) = sum(anomaly_map_rx & normal_map) / sum(normal_map);
    PD(index2) = sum(anomaly_map_rx & anomaly_map) / sum(anomaly_map);
end
end
