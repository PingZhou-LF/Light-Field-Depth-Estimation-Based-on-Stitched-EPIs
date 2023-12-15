%% 直线检测 410s
tic
[Nt,Ns,Ny,Nx,~]=size(LF);
csai=squeeze(LF((Nt+1)/2,(Ns+1)/2,:,:,:));
thresDis=0.1;%EPI区域生长dis0阈值
labelSx={};
for y_ref=1:1:Ny
    for x_ref=1:1:Nx
        [~,~,~,xmin,xmax] =  regionGrowEPI(x_ref,y_ref,LF,thresEPI,dis0,thresDis,'sx');
        xmin=xmin+1;
        xmax=xmax-1;
        if(~isempty(xmin) && ~isempty(xmax)  && xmin<xmax )
            [index_y_temp,index_x_temp]=meshgrid(y_ref,xmin:xmax);
            labelSx{y_ref,x_ref}(:,1)=index_y_temp;
            labelSx{y_ref,x_ref}(:,2)=index_x_temp;
        else
            labelSx{y_ref,x_ref}=[];
        end
    end
end
toc

%得到纹理图
flagEdgeSx=zeros(Ny,Nx);
for i=1:1:Nx
    for j=1:1:Ny
        lenLabel=size(labelSx{j,i},1);
        if(lenLabel<thresLineNum)
            flagEdgeSx(j,i)=1;
        end
    end
end
%  figure;imshow(imoverlay(csai,flagEdgeSx,'b'),[],'InitialMagnification','fit');title('flagEdgeSx')
% 去除线段中的纹理点
for i=1:1:Nx
    for j=1:1:Ny
        temp=labelSx{j,i};
        lenTemp=size(temp,1);%当前线段的元素数
        for k=1:1:lenTemp
            yTemp=labelSx{j,i}(k,1);
            xTemp=labelSx{j,i}(k,2);
            if(flagEdgeSx(yTemp,xTemp)==1)
                [yInd,xInd]=find(temp(:,1)==yTemp & temp(:,2)==xTemp);
                temp(yInd,:)=[];
            end
        end
        labelSx{j,i}=temp;
     end
end
%% 行内merge 生成TLM待选点+SxLine列表 88s
tic
tlmIndSx=cell(Ny,1);
sxLine={};
sxLineInd=1;
for y0=1:1:Ny
    % %行内merge
    linList={};
    listInd=1;
    % figure;imshow(squeeze(LF(5,:,y,:,:)))
    xStart=1;
    while(isempty(labelSx{y0,xStart}) && xStart<Nx)%找到第一个线段长度不为空的像素点
            xStart=xStart+1;
    end
    if(xStart==Nx)
        continue;
    end
    linList{listInd,1}=labelSx{y0,xStart};
    listInd=listInd+1;
    for x=xStart:1:Nx
        %若当前点对应的线段长度为空/小于阈值
        if(isempty(labelSx{y0,x})| size(labelSx{y0,x},1)<lenThres)
            continue;
        else
            %若可以和前一个线段merge,则更新前一个线段
            if(intersect(labelSx{y0,x},linList{listInd-1,1},'rows'))
                linList{listInd-1,1}=union( linList{listInd-1,1},labelSx{y0,x},'rows');
            %若不可以和前一个线段merge,则更新当前线段
            else
                linList{listInd,1}=labelSx{y0,x};
                listInd=listInd+1;
            end
        end
    end
    %linList调出可查看该行的所有直线分别包含哪些pixel
    %保证线段之间无重合
    lenList=size(linList,1);%确定当前行merge后的线段数目
    lenInd=[];
    for i=2:1:lenList
        temp=linList{i,1}(1,2);%当前表头的x坐标
        while(temp<=linList{i-1,1}(end,2) && temp~=0)
            linList{i,1}=linList{i,1}(2:end,:);
            if(isempty(linList{i}))
                temp=0;
            else
                temp=linList{i,1}(1,2);
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
        sxLine{sxLineInd,1}=linList{i,1};
        sxLineInd=sxLineInd+1;
    end
    tlmIndSx{y0}=lenInd;
end
toc
