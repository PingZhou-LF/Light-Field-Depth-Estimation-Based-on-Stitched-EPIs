function x = checkInd(x,Nx)
%CHECKX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
num = length(x);

for i = 1 : num
    if x(i) < 1, x(i)=1; end
    if x(i) > Nx, x(i) = Nx; end
end
end

