function [dpOut1,cwn]= GlobalOptimize1(info)
alpha=info.alpha;
confi= info.c;
cocc=info.cocc;
csai=info.csai;
dis2=info.dis2;
TLM = info.TLM;
dis0 = info.dis0;
disOcc = info.disOcc;
spInfo.cwParams.Th_textureSupress= 0;
%% step1:cw manipulation 5.3.1
cw = confi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tlm
cw(TLM==1)=exp(mean(mean(cw)) - 1);% (5-4)


estVar = stdfilt(dis2(:,:,1),ones(3,3));
estVar= estVar/min(max(estVar(:)),10*median(estVar(:)));
Th_minEstVar=info.Th_minEstVar;
stimes= 10000; % stimes=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5-3
estVarStretch= (estVar<=Th_minEstVar)+ (estVar>Th_minEstVar)*2./(1+exp(stimes*(estVar- Th_minEstVar)));
cwn= cw.*estVarStretch;

%% step2:edge confidence manipulation 5.3.2
occEdgeShrinkage= ones(size(cw,1),size(cw,2));
lcwEdgeShrinkage= ones(size(cw,1),size(cw,2));
tlmEdgeShrinkage= ones(size(cw,1),size(cw,2));
occEdgeEnforce= info.occEdgeEnforce; %5; %5
lcwEdgeEnforce=info.lcwEdgeEnforce;
Th_lcwEdgEnforce= info.Th_lcwEdgEnforce; % 0.1
vals = cw;
mask1 = dis0 ~= disOcc;
occEdgeShrinkage= occEdgeShrinkage.*(mask1).*(1+occEdgeEnforce*cos(pi*vals/2))+... %mask1=1
    occEdgeShrinkage.*(~mask1);% 5-5

lcwEdgeShrinkage= lcwEdgeShrinkage.*(cwn<Th_lcwEdgEnforce).*(1+lcwEdgeEnforce*cos(pi*cwn/2))+...
    lcwEdgeShrinkage.*(cwn>=Th_lcwEdgEnforce);% 5-7

tlmShrinkage= tlmEdgeShrinkage.*(TLM).*2+... %mask1=1
    tlmEdgeShrinkage.*(~TLM);% 5-6

% figure;imshow(lcwEdgeShrinkage,[])
spInfo.occEdge= occEdgeShrinkage.*lcwEdgeShrinkage./tlmShrinkage;% 5-8

%% *******************************
% dis2：待优化的深度图
% cwn：置信度
% csai：中心子孔径图
%  para.alpha=0.0001;%0.0001 alpha相当于(5-1)lamda
% spInfo.occEdge:a(x)
dpOut1 = wls_optimization({dis2}, {cwn} ,...
    csai, alpha, spInfo.occEdge, 0);
% dpOut1(dpOut1>length(spInfo.dpRes))=length(spInfo.dpRes);
% dpOut1(dpOut1<0)=0;
% figure;imshow([gradIn,gradIn.*spInfo.occEdge]);
% mask= sum(TLM,3);
% mask(find(~mask))=nan;
% figure(1);
% subplot(1,2,1);imshow(-info.dis0,[]); title('before')
% subplot(1,2,2);imshow(-dpOut1,[]);title('after')
% [xx,yy] = meshgrid(1:size(dpOut1,1));
% figure(2);surf(xx,yy,dpOut1.*mask);

end

