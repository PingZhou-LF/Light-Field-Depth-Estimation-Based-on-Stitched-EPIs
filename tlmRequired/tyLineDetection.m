%% 直线检测 410s
tic
[Nt,Ns,Ny,Nx,~]=size(LF);
csai=squeeze(LF((Nt+1)/2,(Ns+1)/2,:,:,:));
labelTy={};
for y_ref=1:1:Ny
    for x_ref=1:1:Nx
        [~,~,~,ymin,ymax] =  regionGrowEPI(x_ref,y_ref,LF,thresEPI,dis0,thresDis,'ty');
        ymin=ymin+1;
        ymax=ymax-1;
        if(~isempty(ymin) && ~isempty(ymax) && ymin<ymax )
            [index_y_temp,index_x_temp]=meshgrid(ymin:ymax,x_ref);
            labelTy{y_ref,x_ref}(:,1)=index_y_temp;
            labelTy{y_ref,x_ref}(:,2)=index_x_temp;
        else
            labelTy{y_ref,x_ref}=[];
        end
    end
end
toc

%得到纹理图
flagEdgeTy=zeros(Ny,Nx);
for i=1:1:Nx
    for j=1:1:Ny
        lenLabel=size(labelSx{j,i},1);
        if(lenLabel<thresLineNum)
            flagEdgeTy(j,i)=1;
        end
    end
end
%  figure;imshow(imoverlay(csai,flagEdgeTy,'r'),[],'InitialMagnification','fit');title('flagEdgeTy')
% 去除线段中的纹理点
for i=1:1:Nx
    for j=1:1:Ny
        temp=labelTy{j,i};
        lenTemp=size(temp,1);%当前线段的元素数
        for k=1:1:lenTemp
            yTemp=labelTy{j,i}(k,1);
            xTemp=labelTy{j,i}(k,2);
            if(flagEdgeTy(yTemp,xTemp)==1)
                [yInd,xInd]=find(temp(:,1)==yTemp & temp(:,2)==xTemp);
                temp(yInd,:)=[];
            end
        end
        labelTy{j,i}=temp;
     end
end
%% 行间merge 生成TLM待选点+TyLine列表 88s
tic
tlmIndTy=cell(1,Nx);
tyLine={};
tyLineInd=1;
for x0=1:1:Nx
    % %行内merge
    linList={};
    listInd=1;
    % figure;imshow(squeeze(LF(5,:,y,:,:)))
    yStart=1;
    while(isempty(labelTy{yStart,x0}) && yStart<Ny)%找到第一个线段长度不为空的像素点
            yStart=yStart+1;
    end
    if(yStart==Ny)
        continue;
    end
    linList{listInd,1}=labelTy{yStart,x0};
    listInd=listInd+1;
    for y=yStart:1:Ny
        %若当前点对应的线段长度为空/小于阈值
        if(isempty(labelTy{y,x0})| size(labelTy{y,x0},1)<lenThres)
            continue;
        else
            %若可以和前一个线段merge,则更新前一个线段
            if(intersect(labelTy{y,x0},linList{listInd-1,1},'rows'))
                linList{listInd-1,1}=union( linList{listInd-1,1},labelTy{y,x0},'rows');
            %若不可以和前一个线段merge,则更新当前线段
            else
                linList{listInd,1}=labelTy{y,x0};
                listInd=listInd+1;
            end
        end
    end
    %linList调出可查看该行的所有直线分别包含哪些pixel
    %保证线段之间无重合
    lenList=size(linList,1);%确定当前行merge后的线段数目
    lenInd=[];
    for i=2:1:lenList
        temp=linList{i,1}(1,1);%当前表头的y坐标
        while(temp<=linList{i-1,1}(end,1) && temp~=0)
            linList{i,1}=linList{i,1}(2:end,:);
            if(isempty(linList{i}))
                temp=0;
            else
                temp=linList{i,1}(1,1);
            end
        end
        if temp==0
            break;
        end
    end
    
    %数据整合 
    for i=1:1:lenList
        if(isempty(linList{i}) || size(linList{i},1)<lenThres)
            continue;
        end
        lenInd=[lenInd;linList{i,1}];
        tyLine{tyLineInd,1}=linList{i,1};
        tyLineInd=tyLineInd+1;
    end
    tlmIndTy{x0}=lenInd;
end
toc


