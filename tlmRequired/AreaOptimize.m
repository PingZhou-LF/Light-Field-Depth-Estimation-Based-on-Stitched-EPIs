function [re_depth,tlm] = AreaOptimize( mask,depth )
% 输入：同质区域分类结果，初始深度图
% type=1 线性；type=2 非线性

% 输出：优化后结果
% q:区域个数

[row,col] = size(mask);
re_depth = depth;
%边界优化
[EgPts,mask0] = GetLines(mask,depth);%mask0:边界
if isempty(EgPts)
    re_depth=depth;
    tlm=[];
    return
end
mask0=imfill(mask0,'holes');
tlm = mask0;
nedges = size(EgPts,1);
fsum=@(x)(size(x{:},1));
maxval = arrayfun(fsum,EgPts,'UniformOutput',0);
maxval = max(cell2mat(maxval));
LinesDepth = zeros(maxval,nedges); %每条边界上的深度
LinesCo    = zeros(maxval,2,nedges); %每条边界点的坐标

for ln = 1:nedges
    data = cell2mat(EgPts(ln));
    if(isempty(data))
        re_depth=depth;
        tlm=[];
        return 
    end
    ptnum = size(data,1);
    LinesCo(1:ptnum,1:2,ln) = data(1:ptnum,1:2);
    LinesDepth(1:ptnum,ln) = data(1:ptnum,3);            
    LinesDepth((ptnum+1):end,ln)=inf;
end %赋值

BoardDepth = zeros(row,col);
OptBoardsLines=zeros(row,col);
for t = 1:size(LinesCo,1)
    for p = 1:size(LinesCo,3)
        y = LinesCo(t,1,p);
        x = LinesCo(t,2,p);
        if y==0 || x==0, continue;
        else
            BoardDepth(y,x) = LinesDepth(t,p);
            OptBoardsLines(y,x) = 1;
        end
    end
end %将新的边界值赋给像素点

[xx,yy] = meshgrid(1:col,1:row);
% BoardDepth(~BoardDepth)= NaN;
BoardDepth(BoardDepth==0)= NaN;
if row*col-sum(sum(isnan(BoardDepth)))<=1
    re_depth=depth;
    tlm=[];
    return
end
[fitresult, ~] = createFit(xx, yy, BoardDepth);

for i = 1:col
    for j=1:row
        if mask0(j,i)
           re_depth(j,i) = fitresult.p00+fitresult.p10*i+fitresult.p01*j;
        end
    end
end


a = re_depth.*mask0;
a(a==0)=NaN;
c= depth.*mask0;
c(c==0)=NaN;
[xx,yy] = meshgrid(1:col,1:row);


end 




