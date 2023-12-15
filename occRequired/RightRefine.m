function disRefine=RightRefine(LF,dis_occ,info,maskRight,kRes,sai)
%% 参数输入
% figure;imshow(labeloverlay(gt,3*maskLeft),[],'InitialMagnification','fit');title("maskAll")
[Nt,Ns,~,Nx,~]=size(LF);
Rs=(Ns+1)/2;
dis=dis_occ;
disRefine=dis_occ;
[rightY,rightX]=find(maskRight==1);
tic
for i=1:1:size(rightX)
x=rightX(i);
y=rightY(i);
%% 计算斜率范围
k0 = 1/dis(y,x);
b0 = (Nt*Rs+1)/2-k0*x;%半拼接极图
[kmin,kmax] = calRangeK(k0,b0,Nt*Rs);
%% 拟合最优变化率
costDelta=10000.*ones(201,1);
costInd=1;
for deltaTemp=-0.1:0.001:0.1
    kTemp = k0./(1+(0:1:4)*deltaTemp*k0);%右时(0:1:4)
    xTemp = x+(0:1:4);
    if(xTemp(5)<=Nx && sum(kTemp==0)~=5 &&  sum(isnan(kTemp))==0 )
        [~,~,costTemp1]= indexCalDown(info,kTemp(1),xTemp(1),y,sai(:,Rs:Ns,:));
        [~,~,costTemp2]= indexCalDown(info,kTemp(2),xTemp(2),y,sai(:,Rs:Ns,:));
        [~,~,costTemp3]= indexCalDown(info,kTemp(3),xTemp(3),y,sai(:,Rs:Ns,:));
        [~,~,costTemp4]= indexCalDown(info,kTemp(4),xTemp(4),y,sai(:,Rs:Ns,:));
        [~,~,costTemp5]= indexCalDown(info,kTemp(5),xTemp(5),y,sai(:,Rs:Ns,:));
        costDelta(costInd,1)=costTemp1+costTemp2+costTemp3+costTemp4+costTemp5;
        costInd=costInd+1;
    else
        costDelta(costInd,1)=10000;
        costInd=costInd+1;
    end
end

[~,deltaInd]=min(costDelta);
%figure;plot(costDelta)
delta=-0.1+0.001*(deltaInd-1);
%%
cost=10000.*ones(kRes+1,1);
costInd=1;
kstep = (kmax - kmin) / (kRes - 1);
for indK=0:1:kRes
    k=kmin + indK * kstep;
    kTemp = k./(1+(0:1:4)*delta*k);%flag=2时(0:1:4)
    xTemp = x+(0:1:4);
    if(xTemp(5)<=Nx  && sum(kTemp==0)~=5 &&  sum(isnan(kTemp))==0 )%flag=2时xTemp(5)<=Nx
        [~,~,costTemp1]= indexCalDown(info,kTemp(1),xTemp(1),y,sai(:,Rs:Ns,:));
        [~,~,costTemp2]= indexCalDown(info,kTemp(2),xTemp(2),y,sai(:,Rs:Ns,:));
        [~,~,costTemp3]= indexCalDown(info,kTemp(3),xTemp(3),y,sai(:,Rs:Ns,:));
        [~,~,costTemp4]= indexCalDown(info,kTemp(4),xTemp(4),y,sai(:,Rs:Ns,:));
        [~,~,costTemp5]= indexCalDown(info,kTemp(5),xTemp(5),y,sai(:,Rs:Ns,:));
        cost(costInd,1)=costTemp5*exp(-1) + costTemp4*exp(-9/16) +costTemp3*exp(-4/16) +costTemp2*exp(-1/16) +costTemp1;
        costInd=costInd+1;
    else
        cost(costInd,1)=10000;
        costInd=costInd+1;
    end 
end
[~,costInd]=min(cost);
k=kmin+kstep*(costInd-1);
disTemp=threshold(1./k,info.dis_max,info.dis_min);
disRefine(y,x)=disTemp;
end
toc
% figure;imshow(abs(dis-disRefine),[],'InitialMagnification','fit');
% figure;imshow(dis,[],'InitialMagnification','fit');
% figure;imshow(disRefine,[],'InitialMagnification','fit');