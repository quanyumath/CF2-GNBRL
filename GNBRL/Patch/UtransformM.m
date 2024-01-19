function U = UtransformM(Ycell)
   LL = size(Ycell,1);
   O = cell(LL,1);
   for i = 1:LL
    A = tenmat(Ycell{i},[3]); %square norm
    O{i} = A.data;
    [U{i} D{i} V{i}] = svd(O{i},'econ');
    clear D{i} V{i} A
   end

end