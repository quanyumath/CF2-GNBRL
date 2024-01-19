function [Turecell, Tcell, Maskcell, Scell, sizePatch, index, par, L] = match(Ture, Tensor, Mask, Our_S)

%% parameter setting
par.patsize = 9; %9
par.patnum = 50; % 50
par.step = floor((par.patsize - 1));

sizeData = size(Tensor);
Tpatch = Im2Patch3D(Tensor, par);
sizePatch = size(Tpatch);
[Sel_arr] = nonLocal_arr(sizeData, par); % select seed patchs of groupes
L = length(Sel_arr); % number of groups
Maskpatch = Im2Patch3D(Mask, par);
Turepatch = Im2Patch3D(Ture, par);
Spatch = Im2Patch3D(Our_S, par);
% form a 3rd-order tensor by stacking mode-3 unfolding of  overlapped cubes

%% block matching to find similar cubes
unfoldPatch = Unfold(Tpatch, sizePatch, 2)';
patchXpatch = sum(unfoldPatch.*unfoldPatch, 1);
distenMat = repmat(patchXpatch(Sel_arr), sizePatch(2), 1) + repmat(patchXpatch', 1, L) - 2 * (unfoldPatch') * unfoldPatch(:, Sel_arr);

[~, index] = sort(distenMat);
index = index(1:par.patnum, :);
Tcell = cell(L, 1);
Turecell = cell(L, 1);
Maskcell = cell(L, 1);
Scell = cell(L, 1);
clear unfoldPatch distenMat O patchXpatch

for i = 1:L
    Tcell{i} = Tpatch(:, index(:, i), :);
    Turecell{i} = Turepatch(:, index(:, i), :);
    Maskcell{i} = Maskpatch(:, index(:, i), :);
    Scell{i} = Spatch(:, index(:, i), :);
end