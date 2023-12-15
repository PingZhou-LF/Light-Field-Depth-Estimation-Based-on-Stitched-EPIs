function x = checkInd(x,Nx)
%CHECKX 此处显示有关此函数的摘要
%   此处显示详细说明
num = length(x);

for i = 1 : num
    if x(i) < 1, x(i)=1; end
    if x(i) > Nx, x(i) = Nx; end
end
end

