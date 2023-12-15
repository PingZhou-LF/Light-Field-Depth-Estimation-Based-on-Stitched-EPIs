function [areaGrowTlm,TLM,flagAfterGrow] = areaGrow(mask,img,thres,flag)
%areaGrowTlm:生长后的mask
%TLM:经过填洞、边缘平滑的mask
% 生长clue：
% 生长限制：1、是否归属于其他聚类 2、是否属于图像边界

maxImg=max(max(img));
img=255/maxImg*img;
boundaries=bwboundaries(mask);% figure(1);imshow(mask);
boundY=boundaries{1}(:,1);
boundX=boundaries{1}(:,2);
% [centerY,centerX]=centerSelect(LF,areaAll);
boundNum=size(boundY,1);
areaGrowTlm=mask;
flagBefore=flag;
for i=1:1:boundNum
    xStart=boundX(i);
    yStart=boundY(i);
     [~,label] = regionGrowGray(flag,xStart,yStart,img,thres);
    areaGrowTlm=areaGrowTlm|label;
    flag=flag|label;
%     figure(1);imshow(areaGrowTlm);hold on;plot(xStart,yStart,'r.')
%     title(i)
end


TLM=imfill(areaGrowTlm, 'holes');
kernel = 3;
TLM= medfilt2(TLM,[kernel ,kernel ]);
flagAfterGrow=flag;
% figure;imshow(TLM,[]);