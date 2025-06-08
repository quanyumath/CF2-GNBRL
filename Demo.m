currentFolder = pwd; addpath(genpath(currentFolder));
clear; close all;
clc
load abu-airport-2
mask = map; 

% Airport1--[3,0.3], Airport2--[2,1.2]     100*100*205
% Urban5--[4,0.3]                          100*100*205
% Beach4--[3,1.2]                          150*150*102

% f_show = data(:,:,[37,18,8]);
DataTest = data;
[H, W, Dim] = size(DataTest);
num = H * W;
for i = 1:Dim
    DataTest(:, :, i) = (DataTest(:, :, i) - min(min(DataTest(:, :, i)))) / (max(max(DataTest(:, :, i))-min(min(DataTest(:, :, i)))));
end

numb_dimension = 2;
lambda(1) = 0.01; lambda(2) = 0.01;
X = PCA_img(DataTest, numb_dimension);
A0 = dictionary_A0(X, 'CapLog');

% ----Coarse Stage----
[Our_CA, Our_CL, Our_CS] = acc_GNBRL_solver(X, 50, 0.05, [1, 1, 1]/3, ...
    1.5, lambda, W, A0, 'CapLog');

% ----Fine Stage ----
for k = 1:numb_dimension Mask(:, :, k) = mask; end
[Turecell, Tcell, Mcell, Scell, sizePatch, index, par, LL] = match(X, A0, Mask, Our_CS);
parfor k = 1:LL
    [Our_A{k}, Our_L{k}, Our_S{k}] = acc_GNBRL_solver(Turecell{k}, 50, 0.05, [1, 1, 1]/3, ...
        1.5, lambda, size(Turecell{k}, 2), Tcell{k}, 'CapLog');
end

% Patch Comparison & Replacement
Our_PS = Scell;
for k = 1:LL
    if norm(Our_S{k}(:)-Scell{k}(:)) / norm(Our_S{k}(:)) > 1.2 %0.5,1.2,0.5
        Our_PS{k} = Our_S{k};
    end
end

% reconstruct
Epatch = zeros(sizePatch);
Weight = zeros(sizePatch(1), sizePatch(2));
for i = 1:LL
    Epatch(:, index(:, i), :) = Epatch(:, index(:, i), :) + Our_PS{i};
    Weight(:, index(:, i)) = Weight(:, index(:, i)) + ones(size(Tcell{i}, 1), size(Tcell{i}, 2));
end
[SUNon, ~] = Patch2Im3D(Epatch, Weight, par, size(X)); 
Show_Our = sqrt(sum(SUNon.^2, 3));
Show_Our = (Show_Our - min(Show_Our(:))) / (max(Show_Our(:)) - min(Show_Our(:)));
AUC_Our = AUC(Show_Our, mask); 

