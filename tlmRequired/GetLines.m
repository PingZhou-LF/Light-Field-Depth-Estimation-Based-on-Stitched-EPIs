function [EgPts,region] = GetLines(mask,depth)
%GETLINES 获取mask的各边界直线和点
%mask：=AreaMask
%EgPts

bw = imfill(mask,'holes');
%%%%%%%%%%%%%%%%补边缘
% % 上
ind = find(bw(1,:));
bw(1,min(ind):max(ind))=1;
% 下
ind = find(bw(end,:));
bw(end,min(ind):max(ind))=1;
% 左
ind = find(bw(:,1));
bw(min(ind):max(ind),1)=1;
% 右
ind = find(bw(:,end));
bw(min(ind):max(ind),end)=1;

bw = imfill(bw,'holes');
% figure;imshow(bw)
h = fspecial('average',7);
bw1=imfill(imfilter(double(bwperim(bw)),h)>0,'holes');
% figure;imshow(bw1)
bw2=imerode(bw1,strel('square',9));
% figure;imshow(bw2)

%%%%%%%%%%%%角点检测
C = detectHarrisFeatures(bw2);
C=C.Location;
C=[C(:,2),C(:,1)];%[y,x]

%%%%%%%%边缘交点纳入角点
tmp =zeros(size(bw2));
%上
ind = find(bw2(1,:));
minind = min(ind);
maxind = max(ind);
tmp(1,minind)=1;
tmp(1,maxind)=1;
%下
ind = find(bw2(end,:));
minind = min(ind);
maxind = max(ind);
tmp(end,minind)=1;
tmp(end,maxind)=1;
%左
ind = find(bw2(:,1));
minind = min(ind);
maxind = max(ind);
tmp(minind,1)=1;
tmp(maxind,1)=1;
%右
ind = find(bw2(:,end));
minind = min(ind);
maxind = max(ind);
tmp(minind,end)=1;
tmp(maxind,end)=1;

ind = find(tmp);
[r,c] = ind2sub(size(tmp),ind);
C=[C;[r,c]];


%%%%%%%%%%%%四舍五入哈
newC = C;
for i=1:size(C,1)
    yx = round(C(i,:));
    if bw2(yx(1),yx(2))
        newC(i,:)=round(newC(i,:));
        continue;
    end
    
    nb = [-1,-1;-1,0;-1,1;0,-1;0,0;0,1;1,-1;1,0;1,1];
    nb(:,1)=nb(:,1)+yx(1);
    nb(:,2)=nb(:,2)+yx(2);
%     nb(:,:,1)=[-1,-1,-1;0,0,0;1,1,1]+yx(1);
%     nb(:,:,2)=[-1,0,1;-1,0,1;-1,0,1]+yx(2);
    nb(nb<1)=1;
    nb(nb(:,1)>size(bw2,1))=size(bw2,1);
    nb(nb(:,2)>size(bw2,2))=size(bw2,2);
    d = sqrt((nb(:,1)-C(i,1)).^2+(nb(:,2)-C(i,2)).^2);
    vals = bw2(sub2ind(size(bw2),nb(:,1),nb(:,2)));
    m = double(vals);
    m(~m)=nan;
    d =m.*d;
    [~,ind] = min(d);
    newC(i,:)=nb(ind,:);
end


% imshow(bw2);title('matlab-conner'),
% hold on
% plot(newC(:,2), newC(:,1), 'r*');


%%%%%%%%%%角点顺时针/逆时针排序
tmp=zeros(size(bw2));
tmp(sub2ind(size(bw2),newC(:,1),newC(:,2)))=2;
tmp=tmp+bw2;
[r,c] = ind2sub(size(bw2),find(bw2, 1 ));
B = bwtraceboundary(bw2,[r c],'N');% B 第一维是y 第二维是x
xy = [];%(y,x)
pp=0;
for i=1:size(B,1)
    y = B(i,1); x=B(i,2);
    if tmp(y,x)==2|tmp(y,x)==3 
        xy=[xy;[y,x]];
        tmp(y,x)=4;
        continue;
    end
    
    nb = [-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
    nb(:,1)=nb(:,1)+y;
    nb(:,2)=nb(:,2)+x;
        nb(nb<1)=1;
    nb(nb(:,1)>size(bw2,1),1)=size(bw2,1);
    nb(nb(:,2)>size(bw2,2),2)=size(bw2,1);
 
    
    vals = tmp(sub2ind(size(tmp),nb(:,1),nb(:,2)));
    ind = find(vals==2|vals==3);
    if ~isempty(ind) 
        pp=pp+1;
        xy=[xy;[nb(ind(1),1),nb(ind(1),2)]];
        tmp(nb(ind(1),1),nb(ind(1),2))=4;
    end
end





%%
newxy = xy;
for i=1:size(xy,1)-1
    pt1 = newxy(i,:);
    pt2 = newxy(i+1,:);
    distance = norm(pt1-pt2);
    if distance<=3
        pt = round((pt1+pt2)/2);
        newxy(i,:)=nan;
        newxy(i+1, :)=pt;
    end
end

newxy(isnan(newxy))=[];
if size(newxy,2)~=2
    newxy=reshape(newxy,length(newxy)/2,2);
end
xy=newxy;
%%%%% 
if ~isempty(xy)
    xy=[xy;xy(1,:)];
else
    EgPts=[];region=[];
    return
end
%%%%%%%%%%%
edges = (1:size(xy,1))';
edges = [edges(1:(end-1)),edges(2:end)]; % 起点，终点的像素索引
nedges = size(edges,1);

edgeangles = atan2(xy(edges(:,2),2) - xy(edges(:,1),2), ...
  xy(edges(:,2),1) - xy(edges(:,1),1)); %每条边的角度
k = edgeangles < 0;
edgeangles(k) = edgeangles(k) + 2*pi;
edgeangles = edgeangles/2/pi*360;
nedges = length(edgeangles);
for i = 1:nedges
    edge1 = mod(i-1,nedges)+1;
    edge2 = mod(i,nedges)+1;
    
    angle1 = edgeangles(edge1);
    angle2 = edgeangles(edge2);
    
    if abs(angle1-angle2)<5
        while isnan(edges(edge2,1))
            edge2=edge2+1;
        end
        edges(edge2,1) = edges(edge1,1);
        edges(edge1,:)=NaN;
    end   
end
clear edgeangles
edges(isnan(edges))=[];
% edges=reshape(edges,length(edges)/2,2);
%
if mod(length(edges),2)~=0
    edges=[edges;edges(end,:)];
    edges=reshape(edges,length(edges)/2,[]);
else
    edges=reshape(edges,length(edges)/2,[]);
end

%
nedges = size(edges,1);
%
if isempty(edges)
    EgPts=[];
    region=[];
    return
end

%
edgeangles = atan2(xy(edges(:,2),2) - xy(edges(:,1),2), ...
  xy(edges(:,2),1) - xy(edges(:,1),1)); %每条边的角度
k = edgeangles < 0;
edgeangles(k) = edgeangles(k) + 2*pi;

para1(1) = tan(edgeangles(1));
para2(1)     = xy(edges(1,1),:)*[-para1(end);1];
%  figure(1);imshow(bw2);hold on;
%  plot(xy(edges(1,1),2),xy(edges(1,1),1),'r*'); hold on;
%  hold on;plot([xy(edges(1,1),2),xy(edges(1,2),2)],[xy(edges(1,1),1),xy(edges(1,2),1)],'b-');

ChosenPts=[];
for i=1:nedges
    y1 = xy(edges(i,1),1);
    x1 = xy(edges(i,1),2);
    y2 = xy(edges(i,2),1);
    x2 = xy(edges(i,2),2);
%      hold on;
%      plot(y2,x2,'r*');
%      pause(0.5)
%     distance1 =  abs([x2,y2,1]*[para1(end);-1;para2(end)])/(sqrt((para1(end))^2+1));
%     if distance1<=5, continue; end
     
%     para1 = [para1;tan(edgeangles(i))]; 
%     para2 = [para2;[x1,y1]*[-para1(end);1]];
%      hold on;
%      plot([y1,y2],[x1,x2],'b-');
%      pause(0.5)
    
%     ChosenEgs = [ChosenEgs;edges(i)];
%     ChosenPts= [ChosenPts;xy(edges(i,1),:),xy(edges(i,2),:)];
ChosenPts= [ChosenPts;y1,x1,y2,x2];
end
nedges  = size(ChosenPts,1);

FinalPts = ChosenPts;
% figure(1);imshow(bw2);
for i = 1:nedges
    l1 = mod(i-1,nedges)+1;
    l2 = mod(i,nedges)+1;
    l1y1 = ChosenPts(l1,1); l1x1 = ChosenPts(l1,2);   l1y2 = ChosenPts(l1,3); l1x2 = ChosenPts(l1,4);
    l2y1 = ChosenPts(l2,1); l2x1 = ChosenPts(l2,2);   l2y2 = ChosenPts(l2,3); l2x2 = ChosenPts(l2,4);
    
    if l1x2 == l2x1 && l1y2 ==l2y1, 
        FinalPts(l2,1) = nan; FinalPts(l2,2) = nan;  
        continue; 
    
    end
    
    distance = sqrt((l1x2-l2x1).^2+(l1y2-l2y1).^2);
    if distance<=5
    a = l1y2-l1y1; b=l1x1-l1x2; c = l1x2*l1y1-l1x1*l1y2;
    d = l2y2-l2y1; e=l2x1-l2x2; f = l2x2*l2y1-l2x1*l2y2;
    intery = round((a*f-c*d)/(b*d-a*e));
    interx = round((c*e-b*f)/(b*d-a*e));
    
    FinalPts(l1,3) = intery; FinalPts(l1,4) = interx;
    FinalPts(l2,1) = nan; FinalPts(l2,2) = nan;    
    end
  
%     hold on;
%     plot([l1x1,interx,interx,l2x2],[l1y1,intery,intery,l2y2],'r*-');
    
end
pts =[];
for i=1:2*(size(FinalPts,1))
    r=floor((i-1)/2)+1;
    c=2+(-1).^(mod(i,2));
    yx=FinalPts(r,c:(c+1));
    if sum(isnan(yx)) ~= 0
        continue;
    end
    pts=[pts;yx];
end
if isempty(pts)
    EgPts=[];
    region=[];
    return
end
pts=[pts;pts(1,:)];

% figure; imshow(labeloverlay(CSAI,bw2));
% for i = 1:size(pts,1)-1
%     hold on; 
%     plot([pts(i:i+1,2)],[pts(i:i+1,1)],'r*-');
% end

nedges = size(pts,1)-1;
EgPts = cell(nedges,1);
region = zeros(size(bw2));
for i = 1:nedges
    [tmpy,tmpx]=bresenham(pts(i,1),pts(i,2),pts(i+1,1),pts(i+1,2));
    tmp = zeros(length(tmpx),3);
    v1=tmpy;v2=tmpx;
    siz=size(mask);
    if( min(v1(:)) < 1 || max(v1(:)) > siz(1) || ...
        min(v2(:)) < 1 || max(v2(:)) > siz(2) )
        continue;
    end
    region(sub2ind(size(mask),tmpy,tmpx))=1;
    tmp(:,1) = tmpy;
    tmp(:,2) = tmpx;
    tmp(:,3) = depth(sub2ind(size(depth),tmpy,tmpx));
    EgPts(i) = mat2cell(tmp,size(tmp,1),size(tmp,2));
end


end
