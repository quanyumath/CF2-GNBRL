function Boxplot_AG(mask, Methods_E, Xtick_label)

% The box plot can intuitively reflect the anomaly–background separability.
% The bigger the gap between the abnormal box and the background box,
% the higher the anomaly–background separation will be. The shorter the
% background box, the more background information can be effectively suppressed.

% 说明：
% mask--------真实的异常矩阵
% Methods_E---各种方法得到的异常矩阵
% xlab--------各种方法的名称

label_value = mask(:);
ind_tar = find(label_value == 1);
ind_bac = find(label_value == 0);
num_targ = length(ind_tar);
num_back = length(ind_bac);


X_targ = []; X_back = [];
g1_targ = []; g1_back = [];

num_meth = length(Methods_E);
for k = 1:num_meth
    R = Methods_E{k}; R = normalize(R(:), 'range');
    targ = R(ind_tar)'; back = R(ind_bac)';
    X_targ = [X_targ; targ]; X_back = [X_back; back];
    g1_targ = [g1_targ; k * ones(1, num_targ)];
    g1_back = [g1_back; k * ones(1, num_back)];
end
X_targ = X_targ'; X_back = X_back';
X = [X_targ(:); X_back(:)]; X = X(:);
g1_targ = g1_targ'; g1_back = g1_back';
g1 = [g1_targ(:); g1_back(:)]; g1 = g1(:);
g2 = [ones(num_meth*num_targ, 1); 2 * ones(num_meth*num_back, 1)]; g2 = g2(:);
positions = [[1:num_meth], [1:num_meth] + 0.3];

%%
figure(2);
% 相关属性['plotstyle','compact']['colorgroup',g2,]['color','rk']
bh = boxplot(X, {g2, g1}, 'whisker', 10000, 'colorgroup', g2, 'symbol', '.', 'outliersize', 4, 'widths', 0.2, 'positions', positions, 'color', 'rb');
%grid on; 添加网格线
set(bh, 'LineWidth', 1.5)
ylabel('Detection test statistic range');

% grid on
% set(gca,'YLim',[0,0.5],'gridLineStyle', '-.');

% ylim([0,0.0065])  % 用于y轴的坐标轴显示范围的控制

Xtick_pos = [1:num_meth] + 0.15; % 确定label显示的位置
set(gca, 'XTickLabel', Xtick_label, 'XTick', Xtick_pos); %可以设置字体属性['fontsize',15]
xtickangle(15)% 旋转标签角度

%%
h = findobj(gca, 'Tag', 'Outliers');
delete(h)
%legend(findobj(gca, 'Tag', 'Box'), 'Background', 'Anomaly', 'Orientation', 'horizon')

%% 最大值与最小值（箱须至高与至低点：whisker 为0-100%）
p_targ = prctile(X_targ, [0, 100]);
p_back = prctile(X_back, [0, 100]);
% p_targ = prctile(X_targ,[10 90]);
% p_back = prctile(X_back,[10 90]);
p = [];
for i = 1:num_meth
    p = [p, p_targ(:, i), p_back(:, i)];
end

% 箱子的上边缘与下边缘 (异常、背景区域10% 与 90% 统计)
q_targ = quantile(X_targ, [0.1, 0.9]);
q_back = quantile(X_back, [0.1, 0.9]);
% q_targ = quantile(X_targ,[0.09 0.81]);
% q_back = quantile(X_back,[0.09 0.81]);
q = [];
for i = 1:num_meth
    q = [q, q_targ(:, i), q_back(:, i)];
end

h = flipud(findobj(gca, 'Tag', 'Box'));
for j = 1:length(h)
    q10 = q(1, j);
    q90 = q(2, j);
    set(h(j), 'YData', [q10, q90, q90, q10, q10]);
end

% Replace upper end y value of whisker
h = flipud(findobj(gca, 'Tag', 'Upper Whisker'));
for j = 1:length(h)
    %     ydata = get(h(j),'YData');
    %     ydata(2) = p(2,j);
    %     set(h(j),'YData',ydata);
    set(h(j), 'YData', [q(2, j), p(2, j)]);
end

% Replace all y values of adjacent value
h = flipud(findobj(gca, 'Tag', 'Upper Adjacent Value'));
for j = 1:length(h)
    %     ydata = get(h(j),'YData');
    %     ydata(:) = p(2,j);
    set(h(j), 'YData', [p(2, j), p(2, j)]);
end

% Replace lower end y value of whisker
h = flipud(findobj(gca, 'Tag', 'Lower Whisker'));
for j = 1:length(h)
    %     ydata = get(h(j),'YData');
    %     ydata(1) = p(1,j);
    set(h(j), 'YData', [q(1, j), p(1, j)]);
end

% Replace all y values of adjacent value
h = flipud(findobj(gca, 'Tag', 'Lower Adjacent Value'));
for j = 1:length(h)
    %     ydata = get(h(j),'YData');
    %     ydata(:) = p(1,j);
    set(h(j), 'YData', [p(1, j), p(1, j)]);
end

end