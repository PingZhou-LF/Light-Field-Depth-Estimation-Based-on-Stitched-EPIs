function x = threshold(x,maxX,minX)
% 数据范围限定
for i=1:1:size(x,1)
if x(i,1) < minX, x(i,1)=minX; end
if x(i,1) > maxX, x(i,1) =maxX; end
end
