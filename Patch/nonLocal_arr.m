function  [SelfIndex_arr]  =  nonLocal_arr(sizeD, par)
% -SelfIndex_arr is the index of keypatches in the total patch index array
TempR    =sizeD(1)-par.patsize+1;
TempC    =sizeD(2)-par.patsize+1;
R_GridIdx=1:par.step:TempR;
R_GridIdx=[R_GridIdx R_GridIdx(end)+1:TempR];
C_GridIdx=1:par.step:TempC;
C_GridIdx=[C_GridIdx C_GridIdx(end)+1:TempC];

temp=1:TempR*TempC;
temp=reshape(temp,TempR,TempC);
SelfIndex_arr = temp(R_GridIdx,C_GridIdx);
SelfIndex_arr = SelfIndex_arr(:)';
end

