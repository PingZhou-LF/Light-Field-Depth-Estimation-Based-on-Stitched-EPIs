function [l,tlmMask,tlmAll]=tlmSeg(csai,TLM,hs,hr,ratioThres)
%l:每个像素所在的类标签
%tlmMask：满足区域大小条件的TLM区域
%tlmAll：所有tlm的mask（便于作图
[Ny,Nx,~]=size(csai);
imgTLM=zeros(Ny,Nx,3);
for i=1:1:3
    temp=csai(:,:,i);
    temp(TLM==0)=0;
    imgTLM(:,:,i)=temp;
end
% 提取mask图像
V=rgb2gray(csai);
maxV=max(max(V));
V=255/maxV*V;
% meanshift聚类% l=meanshift(V,dis0,hs,hr,hr) 
l = meanshsegm(V,hs,hr) ;
l(TLM==0)=0;

numAll=sum(sum((TLM==1)));
sort_initial=tabulate(l(:));%value count percent
sort_initial=sort_initial(2:end,:);%去掉标号为0的(非TLM区域)
sort_result=sortrows(sort_initial,2);%按照count比较数值并按升序排序
sort_result(:,2)=sort_result(:,2)/numAll;%计算每个label对应的
index=sort_result(:,2)>ratioThres;%找到大比例区域
labelArea=sort_result(index,1);%找到大比例区域的label

tlmMask=cell(size(labelArea,1),1);
tlmAll=zeros(Ny,Nx);
[Ny,Nx]=size(l);
for i=1:1:size(labelArea,1)
    maskTemp=zeros(Ny,Nx);
    maskTemp(l==labelArea(i))=1;
    tlmMask{i}=maskTemp;
    tlmAll(l==labelArea(i))=1;
end