function diffH = diffH3(X, size)

n1 = size(1); n2 = size(2); n3 = size(3);

dfx = diff(X{1}, 1, 1); dfy = diff(X{2}, 1, 2); dfz = diff(X{3}, 1, 3);

dfxT = zeros(size); dfyT = zeros(size); dfzT = zeros(size);
dfxT(2:end,:,:) = -dfx; dfxT(1,:,:) = X{1}(end,:,:)-X{1}(1,:,:); 
dfyT(:,2:end,:) = -dfy; dfyT(:,1,:) = X{2}(:,end,:)-X{1}(:,1,:); 
dfzT(:,:,2:end) = -dfz; dfzT(:,:,1) = X{3}(:,:,end)-X{3}(:,:,1); 

diffH = dfxT+dfyT+dfzT;
end