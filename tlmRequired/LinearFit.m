function [re_dis]=LinearFit(csai,tlm,dis0,areaThres,musterThres,hs2,hr2,ratioThres)
%tlm:768*768单mask

[Ny,Nx,~]=size(csai);
%孔洞填充
tlmFill=imfill(tlm, 'holes');
% figure;imshow(tlmFill);
%边缘平滑 去除毛刺
kernel = 3;
tlmSmooth= medfilt2(tlmFill,[kernel ,kernel ]);
% figure;imshow(tlmSmooth);
%判断是否连通【有的tlmMask聚类后不再是连通域
[tlmLabel,tlmNum]=bwlabel(tlmSmooth,4);
% figure;imshow(labeloverlay(csai,tlmLabel),[],'InitialMagnification','fit');
disTemp=dis0;
if tlmNum>1
%多连通区域的聚类和优化
    for j=1:1:tlmNum
        if(size(find(tlmLabel==j),1)>areaThres)%表示当前连通域满足区域判定阈值
            mask=tlmLabel;mask(tlmLabel~=j)=0;%提取当前的mask
            if(size(find(tlmLabel==j),1)>musterThres)%表示当前连通域需要进行聚类分析
                r=csai(:,:,1);g=csai(:,:,2);b=csai(:,:,3);
                r(tlmLabel~=j)=0;g(tlmLabel~=j)=0;b(tlmLabel~=j)=0;
                rgb=zeros(Ny,Nx,3);rgb(:,:,1)=r;rgb(:,:,2)=g;rgb(:,:,3)=b;
                hsv = rgb2hsv(rgb);V = hsv(:,:,3); maxV=max(max(V));V=255/maxV*V;
                %tlmClusterLabel = meanshsegm(V,hs2,hr2) ;
                 l = meanshsegm(V,hs2,hr2) ;%对当前连通区域进行聚类
                 numAll=sum(sum((tlmLabel==j)));
                 sort_initial=tabulate(l(:));%value count percent
                 sort_initial=sort_initial(2:end,:);%去掉标号为0的(非TLM区域)
                 sort_result=sortrows(sort_initial,2);%按照count比较数值并按升序排序
                 sort_result(:,2)=sort_result(:,2)/numAll;%计算每个label对应的
                 index=sort_result(:,2)>ratioThres;%找到大比例区域
                 labelArea=sort_result(index,1);%找到大比例区域的label
                 for i=1:1:size(labelArea,1)%对每个聚类做深度优化
                    maskCur=mask;
                    maskCur(l~=labelArea(i))=0;
                    %聚类后处理
                    maskCur=imfill(maskCur, 'holes');%填充空洞
                    maskCur= medfilt2(maskCur,[kernel ,kernel ]);%去毛刺
                    if size(find(maskCur==1),1)>areaThres
                        [re_depth,~] = AreaOptimize( maskCur,disTemp,csai );%深度优化
                        disTemp=re_depth;%深度更新
                    end
                 end
            else%表示当前连通域可直接进行深度优化
                [re_depth,~] = AreaOptimize( mask,disTemp,csai );%深度优化
                disTemp=re_depth;%深度更新
            end
        end
    end
else
    
%单连通域的聚类和优化
if(size(find(tlmSmooth),1)>areaThres)%表示当前连通域满足区域判定阈值
    mask=tlmSmooth;%提取当前的mask
    if(size(find(tlmSmooth),1)>musterThres)%表示当前连通域需要进行聚类分析 
        r=csai(:,:,1);g=csai(:,:,2);b=csai(:,:,3);
        r(tlmSmooth==0)=0;g(tlmSmooth==0)=0;b(tlmSmooth==0)=0;
        rgb=zeros(Ny,Nx,3);rgb(:,:,1)=r;rgb(:,:,2)=g;rgb(:,:,3)=b;
        hsv = rgb2hsv(rgb);V = hsv(:,:,3); maxV=max(max(V));V=255/maxV*V;
        if sum(sum(isnan(V)))>0
            re_dis=dis0;
            return
        end
        l = meanshsegm(V,hs2,hr2) ;
        l(mask==0)=0;
        numAll=sum(sum((mask==1)));
        sort_initial=tabulate(l(:));%value count percent
        sort_initial=sort_initial(2:end,:);%去掉标号为0的(非TLM区域)
        sort_result=sortrows(sort_initial,2);%按照count比较数值并按升序排序
        sort_result(:,2)=sort_result(:,2)/numAll;%计算每个label对应的
        index=sort_result(:,2)>ratioThres;%找到大比例区域
        labelArea=sort_result(index,1);%找到大比例区域的label
        for i=1:1:size(labelArea,1)%对每个聚类做深度优化
            maskCur=mask;
            maskCur(l~=labelArea(i))=0;
            %聚类后处理
            maskCur=imfill(maskCur, 'holes');%填充空洞
            maskCur= medfilt2(maskCur,[kernel ,kernel ]);%去毛刺
            if size(find(maskCur==1),1)>areaThres
            [re_depth,~] = AreaOptimize( maskCur,disTemp,csai );%深度优化
            disTemp=re_depth;%深度更新
            end
        end
    else
        [re_depth,~] = AreaOptimize( mask,disTemp,csai );%深度优化
        disTemp=re_depth;%深度更新
    end
end
end

re_dis=disTemp;
