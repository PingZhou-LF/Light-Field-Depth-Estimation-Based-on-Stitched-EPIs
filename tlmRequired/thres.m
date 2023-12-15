function xNew=thres(x,y)
% 区间到区间映射
%x:原数据
%y：映射后的数据区间
xmax=max(max(x));
xmin=min(min(x));
ymax=max(max(y));
ymin=min(min(y));
k=(ymax-ymin)/(xmax-xmin);
b=ymax-k*xmax;
xNew=k.*x+b;