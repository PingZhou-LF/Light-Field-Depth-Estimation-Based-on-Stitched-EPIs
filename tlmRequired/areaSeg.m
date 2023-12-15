%% 图像边缘信息提取
[Ny,Nx,~]=size(csai);
csaiHSV=rgb2hsv(csai);
S=csaiHSV(:,:,2);
sEdge=zeros(Ny,Nx);
if type==1 %边界明显
    sEdge(abs(gradient(S))>0.05)=1;%S图边界
else%边界不明显
    sEdge(abs(gradient(histeq(S)))>0.05)=1;%S图边界
end
B=[1 1 1 
0 1 0
1 1 1];
sEdgeDilate=imdilate(sEdge,B);%S图边界膨胀
%% 区域提取
areaThres=5*5;
imgTlm=imgMerge&~sEdgeDilate;
tlmFill=imfill(imgTlm,'holes');
tlmMaskAll = bwareaopen(tlmFill,areaThres,4);%去除小面积连通区域
[tlmLabel,bwNum]=bwlabel(tlmMaskAll,4);
% 
% figure;
% subplot(2,2,1);imshow(sEdge);title('S图边界检测')
% subplot(2,2,2);imshow(sEdgeDilate);title('S图边界填洞+膨胀')
% subplot(2,2,3);imshow(imgTlm);title('TLM-S图边界修正')
% subplot(2,2,4);imshow(tlmMaskAll);title('TLM-Final')
% figure;
% subplot(1,2,1);imshow(labeloverlay(csai,tlmLabel),[],'InitialMagnification','fit');title('连通域提取')
% subplot(1,2,2);imshow(tlmLabel,[],'InitialMagnification','fit');title('连通域提取')

%% 区域分割

hs=30;%聚类空间阈值30
hr=1;%聚类色彩阈值4
areaThres=10*10;%判断是否是TLM 10*10
clusterThres=100*100;%判断是否是需要聚类的TLM100*100
threColor=4;%区域生长阈值4

flagAll=zeros(Ny,Nx);%用于区域生长
tlmMaskFinal=zeros(Ny,Nx);%保存label
IndMaskFinal=1;
%flag初始化（+边缘
if type==1 %边界明显
    flagAll(abs(gradient(S))>0.05)=1;%S图边界
else%边界不明显
    flagAll(abs(gradient(histeq(S)))>0.05)=1;%S图边界
end

flagTemp=zeros(bwNum,1);

for i= 1:1:bwNum
    %判断是否属于待优化区域
    if(sum(sum(tlmLabel==i))>clusterThres)%需要进行聚类
        %% 聚类
        if(type==3)
            V=histeq(csaiHSV(:,:,3));
        else
            V=csaiHSV(:,:,3);
        end
        V(tlmLabel~=i)=0;
        maxImg=max(max(V));
        V=255/maxImg*V;
        labelTemp= meanshsegm(V,hs,hr) ;
        % figure;imshow(labelTemp,[],'InitialMagnification','fit') ;
        %% 聚类区域筛选
        sort_initial=tabulate(labelTemp(:));%value count percent
%         sort_initial=sort_initial(2:end,:);%去掉标号为1的(背景)
        sort_result=sortrows(sort_initial,2);%按照count比较数值并按升序排序
        sort_result=sort_result(1:end-1,:);%去掉标号为1的(背景)
        index=sort_result(:,2)>areaThres;
        labelArea=sort_result(index,1);%判断是否满足区域条件
        maskTemp=zeros(Ny,Nx,size(labelArea,1));
        for j=1:1:size(labelArea,1)
            mask_=zeros(Ny,Nx);
            mask_(labelTemp==labelArea(j))=1;
            mask_=imfill(mask_,'holes');
            maskTemp(:,:,j)=mask_;
%         figure(1);imshow(labeloverlay(csai,maskTemp(:,:,j)),[]);pause(1)
        end
        %% 聚类区域生长

        maskGrow=zeros(Ny,Nx,size(labelArea,1));
        maskAllGrow=zeros(Ny,Nx);%便于作图

        for j=1:1:size(labelArea,1)
            mask_= maskTemp(:,:,j);%取mask
            if(sum(V(mask_==1))==0)
                continue;%去除背景
            end
            
            [~,TLM,flagAfterGrow] = areaGrow(mask_,csaiHSV(:,:,3),threColor,flagAll);%区域生长
            flagAll=flagAfterGrow;
            maskGrow(:,:,j)=TLM;
            maskAllGrow(TLM==1)=1;        
            tlmMaskFinal(TLM==1)=IndMaskFinal;
            IndMaskFinal=IndMaskFinal+1;
        end


    else%不需要进行聚类
        if(sum(sum(tlmLabel==i))>areaThres)
            maskTemp=tlmMaskAll;
            maskTemp(tlmLabel~=i)=0;
            maskTemp=imfill(maskTemp, 'holes');
            kernel = 3;
            maskTemp= medfilt2(maskTemp,[kernel ,kernel ]);
            flagAll(maskTemp==1)=1;%更新flagAll
            tlmMaskFinal(maskTemp==1)=IndMaskFinal;
            IndMaskFinal=IndMaskFinal+1;
        end
    end
    X=[num2str(i),'/',num2str(bwNum)];
end

% figure;
% subplot(1,2,1);imshow(labeloverlay(csai,tlmMaskFinal),[],'InitialMagnification','fit');title('TLM检测结果')
% subplot(1,2,2);imshow(tlmMaskFinal,[],'InitialMagnification','fit');title([num2str(hs),',',num2str(hr),',',num2str(threColor)]);
%% 区域修正
se = strel('square',3);
%为了解决聚类所造成的不连通现象
[labelTLM,~] = unique(tlmMaskFinal);
labelTLM=labelTLM(2:end);%去除label=0（非同质区域）
tlmRefine=zeros(Ny,Nx);
IndRefine=1;
for i=1:1:max(labelTLM)
    maskTemp=zeros(Ny,Nx);
    maskTemp(tlmMaskFinal==i)=1;
    %去除毛刺、小连接
    maskOpen=imopen(maskTemp,se);
    maskEx = bwareaopen(maskOpen,areaThres,8) ;
%     figure;
%     subplot(1,3,1);imshow(maskTemp);title('原mask')
%     subplot(1,3,2);imshow(maskOpen);title('开操作')
%     subplot(1,3,3);imshow(maskEx);title('去除小连通域')
    [lTemp,numTemp]=bwlabel(maskEx,8);%figure;imshow(lTemp,[]);
    for j=1:1:max(max(lTemp))
        tlmRefine(lTemp==j)=IndRefine;
        IndRefine=IndRefine+1;
    end   
end
% figure;
% subplot(1,2,1);imshow(labeloverlay(csai,tlmRefine),[],'InitialMagnification','fit');
% subplot(1,2,2);imshow(tlmRefine,[],'InitialMagnification','fit');title('TLM优化结果')


paraTLM.type=type;
paraTLM.B=B;
paraTLM.se=se;
paraTLM.hs=hs;%聚类空间阈值
paraTLM.hr=hr;%聚类色彩阈值
paraTLM.areaThres=areaThres;%判断是否是TLM
paraTLM.clusterThres=clusterThres;%判断是否是需要聚类的TLM
paraTLM.threColor=threColor;%区域生长阈值


