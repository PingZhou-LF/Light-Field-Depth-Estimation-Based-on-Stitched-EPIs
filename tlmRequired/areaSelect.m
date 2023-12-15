function [areaAll,areaMask] = areaSelect(LF,areaThres,imgSx,imgTy,imgMerge)
[Nt,Ns,Ny,Nx,~]=size(LF);
csai=squeeze(LF((Nt+1)/2,(Ns+1)/2,:,:,:));
contourSx = bwperim(imgSx);
contourTy = bwperim(imgTy);
I_gray = rgb2gray(csai);
contourCsai= edge(I_gray, 'log');
% figure(5);
% subplot(1,3,1);imshow(contourSx,[]);title('sx边缘检测')
% subplot(1,3,2);imshow(contourTy,[]);title('ty边缘检测')
% subplot(1,3,3);imshow(contourCsai,[]);title('csai边缘检测')

imFinal=imgMerge &(~contourSx)&(~contourCsai);
% imFinal=imgMerge &(~contourTy)&(~contourCsai);
figure;
subplot(1,3,1);imshow(imFinal,[]);title('imFinal')
%筛选连通区域
[bwLabel,bwNum]=bwlabel(imFinal,4);
bwAreaNum=zeros(bwNum,1);%每个连通区域的像素数目统计
for ii=1:bwNum
    bwAreaNum(ii)=sum(sum(bwLabel==ii));
end
ind=find(bwAreaNum>areaThres);
indNum=size(ind,1);
areaAll=zeros(Ny,Nx);
areaMask=zeros(Ny,Nx,indNum);
mask=zeros(Ny,Nx,indNum);
for i=1:1:indNum
    areaAll=(imfill(bwLabel==ind(i), 'holes')|areaAll);
    %染色
    areaMask(:,:,i)=imfill(bwLabel==ind(i), 'holes'); 
    mask(:,:,i) = areaMask(:,:,i)*i;
end
subplot(1,3,2);imshow(labeloverlay(csai,sum(mask,3)),[],'InitialMagnification','fit')
subplot(1,3,3);imshow(areaAll,[],'InitialMagnification','fit');title('areaAll')

