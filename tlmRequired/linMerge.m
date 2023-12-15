function [imgMerge,imgSx,imgTy] = linMerge(LF,sxLine,tyLine)
%
[~,~,Ny,Nx,~]=size(LF);
sxImg=zeros(Ny,Nx);
tyImg=zeros(Ny,Nx);
sxIndImg=zeros(Ny,Nx);
tyIndImg=zeros(Ny,Nx);
for i=1:1:size(sxLine,1)
    tempInd=sxLine{i,1};
    sxImg(tempInd(:,1),tempInd(:,2))=1;%01
    sxIndImg(tempInd(:,1),tempInd(:,2))=i;%编号
end
% figure(1);
% subplot(1,2,1);imshow(sxImg,[]);title('HTLL')
% subplot(1,2,2);imshow(sxIndImg,[]);title('sx直线编号')

for i=1:1:size(tyLine,1)
    tempInd=tyLine{i,1};
    tyImg(tempInd(:,1),tempInd(:,2))=1;
    tyIndImg(tempInd(:,1),tempInd(:,2))=i;
end
% figure(2);
% subplot(1,2,1);imshow(tyImg,[]);title('VTLL')
% subplot(1,2,2);imshow(tyIndImg,[]);title('ty直线编号')


sxIndImgRe=sxIndImg;
tyIndImgRe=tyIndImg;
imgCross = zeros(Ny,Nx);%交点01
imgCross(sxImg & tyImg ==1)=1;
sxIndImgRe(imgCross~=1)=0;%sx交点编号
tyIndImgRe(imgCross~=1)=0;%ty交点编号
[Csx,~,~] = unique(sxIndImgRe); 
imgSx=zeros(Ny,Nx);
for i=2:1:size(Csx,1)
    imgSx(sxIndImg==Csx(i))=1;%交点扩张
end
% figure(3);imshow(imgSx,[],'InitialMagnification','fit');title('sx交点扩张')

[Cty,~,~] = unique(tyIndImgRe); 
imgTy=zeros(Ny,Nx);
for i=2:1:size(Cty,1)
    imgTy(tyIndImg==Cty(i))=1;
end
% figure(3);imshow(imgTy,[],'InitialMagnification','fit');title('ty交点扩张')
imgMerge=zeros(Ny,Nx);
imgMerge(imgTy|imgSx==1)=1;
% figure(4);imshow(imgMerge,[],'InitialMagnification','fit'); title('Merge')