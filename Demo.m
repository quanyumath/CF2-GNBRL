currentFolder = pwd; addpath(genpath(currentFolder));

clear; close all;
clc

jj = 1;
if jj == 1
    load abu-airport-1
    numb_dimension = 3; u = 0.65;
elseif jj == 2
    load abu-airport-2
    numb_dimension = 2; u = 0.78;
elseif jj == 3
    load abu-urban-5
    numb_dimension = 4; u = 0.65;
elseif jj == 4
    load abu-beach-4
    numb_dimension = 3; u = 1.2;
end

mask = map;

f_show = data(:, :, [37, 18, 8]);
DataTest = data;
[H, W, Dim] = size(DataTest);
num = H * W;
for i = 1:Dim
    DataTest(:, :, i) = (DataTest(:, :, i) - min(min(DataTest(:, :, i)))) / (max(max(DataTest(:, :, i))-min(min(DataTest(:, :, i)))));
end
Y = reshape(DataTest, num, Dim)';

X = PCA_img(DataTest, numb_dimension);
A0 = dictionary_A0(X, 'CapLog');
% imshow(X(:,:,3))

%% GNBRL
opts = [];
opts.maxit = 50;
opts.tol = 1e-2;
opts.A = A0;
opts.psi = 'CapL1';
tic
[Our_CA, Our_CL, Our_CS] = acc_GNBRL_solver(X, opts);
%time_Our = toc
AUC_Our = AUC(Our_CS, mask) %
Show_Our = sqrt(sum(Our_CS.^2, 3));
Show_Our = (Show_Our - min(Show_Our(:))) / (max(Show_Our(:)) - min(Show_Our(:)));


% ----Fine Stage ----
for k = 1:numb_dimension Mask(:, :, k) = mask; end
[Turecell, Tcell, Mcell, Scell, sizePatch, index, par, LL] = match(X, A0, Mask, Our_CS);
parfor k = 1:LL
    [Our_A{k}, Our_L{k}, Our_S{k}] = CF2_GNBRL_solver(Turecell{k}, Tcell{k});
end
% Patch Comparison & Replacement
Our_PS = Scell;
for k = 1:LL
    if norm(Our_S{k}(:)-Scell{k}(:)) / norm(Our_S{k}(:)) > u %0.5,1.2,0.5
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
[SUNon, ~] = Patch2Im3D(Epatch, Weight, par, size(X)); % imshow(SUNon(:,:,2))
time_Our = toc
Show_CF_Our = sqrt(sum(SUNon.^2, 3));
Show_CF_Our = (Show_CF_Our - min(Show_CF_Our(:))) / (max(Show_CF_Our(:)) - min(Show_CF_Our(:)));
AUC(SUNon, mask) % imshow(Show_Our) OL = Our_S(:,:,60); imshow(OL/max(OL(:))+0.5)         imshow(mask)

%% ROC
Re_tensor = {};
Re_tensor{1} = Show_Our; Re_tensor{2} = Show_CF_Our;
methodname = {'GNBRL', 'CF2-GNBRL'};
clear Methods_E
clear xlab
for i = 1:2
    Methods_E{i} = Re_tensor{i};
    xlab{i} = methodname{i};
end

%% ROC
ROC(Methods_E, mask, xlab)

%% Boxplot
Boxplot_AG(mask, Methods_E, methodname)
