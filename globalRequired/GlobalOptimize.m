function [dpOut1,cwn]= GlobalOptimize(para)
alpha=para.alpha;
c = para.c;
cocc=para.cocc;
csai=para.csai;
dis2=para.dis2;
TLM = para.TLM;
spInfo.cwParams.Th_textureSupress= 0;


cw = c.*cocc;
% cw = confi;
cwTLM = zeros(size(cw));
% centralView= im2double(spInfo.centralView);

% 1. low confidence occlusion border shrinkage
for i=1:size(TLM,3)
    TLM(:,:,i) = imerode( TLM(:,:,i),strel('square',5));
end


for i=1:size(TLM,3)
    mask = TLM(:,:,i);
    boundary = edge(TLM(:,:,i),'canny');
    ind = find(boundary);
    meanconfi = mean(cw(ind));
    if meanconfi < 1
        meanconfi = 1;
    end
    cwboundary = boundary.*cw;
    cwboundary(cwboundary<meanconfi) = meanconfi;
    cwTLM = cwTLM+cwboundary.*boundary+meanconfi*(mask-boundary);
end
% cwn = cw+cwTLM;

% % 2. high confidece intricate structure preservation 保留高置信度区域
% hcMask= spInfo.dpWeightSPcor(:,:,1)>0.5; 
% x= (~imopen(hcMask,strel('square',5))).*hcMask; % figure;imshow(x);
% dpStretch= (~x).*dpStretch+ x; % figure;imshow(dpStretch);

% % 3. high est. variation region supression eq10
estVar = stdfilt(dis2(:,:,1),ones(3,3));
estVar= estVar/min(max(estVar(:)),10*median(estVar(:)));
Th_minEstVar=para.Th_minEstVar;
stimes= 10000; % stimes=10000;
% mask1 = estVar<=Th_minEstVar;
% mask1=imerode(mask1,strel('square',5));
% mask1=imdilate(mask1,strel('square',5));
estVarStretch= (estVar<=Th_minEstVar)+ (estVar>Th_minEstVar)*2./(1+exp(stimes*(estVar- Th_minEstVar)));
% estVarStretch= mask1+ (~mask1)*2./(1+exp(stimes*(estVar- Th_minEstVar))); 
% figure;imshow(estVarStretch);

cwn= cw.*estVarStretch;

mask = sum(TLM,3);
cwn(find(mask))=0.8;

%% edge confidence manipulation 
occEdgeShrinkage= ones(size(cw,1),size(cw,2));
% % % Th_minResp= 1.5;
% % % occEdgeShrinkage= (spInfo.minResp<Th_minResp)*2./(1+exp(-1*(spInfo.minResp-Th_minResp)))+ (spInfo.minResp>=Th_minResp);
occEdgeEnforce= para.occEdgeEnforce; %5; %5
lcwEdgeEnforce= para.lcwEdgeEnforce; % 2;
 % 2;

mask1 = cocc<1;
vals = cocc;
mask1=imdilate(mask1,strel('square',5));
vals(find(cocc>=1))= 0;

occEdgeShrinkage= occEdgeShrinkage.*(mask1).*(1+occEdgeEnforce*cos(pi*vals/2))+...
                    occEdgeShrinkage.*(~mask1);
%       

% Th_lcwEdgEnforce= 0.2; % 0.1
% occEdgeShrinkage= occEdgeShrinkage.*(cwn<Th_lcwEdgEnforce).*(1+lcwEdgeEnforce*cos(pi*cwn/2))+...
%                     occEdgeShrinkage.*(cwn>=Th_lcwEdgEnforce);
spInfo.occEdge= occEdgeShrinkage; 
% figure; imshow(occEdgeShrinkage);
% figure; imshow(spInfo.occEdge/5);

dpOut1 = wls_optimization({dis2}, {cwn} ,...
        csai, alpha, spInfo.occEdge, 0);
% dpOut1(dpOut1>length(spInfo.dpRes))=length(spInfo.dpRes);
% dpOut1(dpOut1<0)=0;
% figure;imshow([gradIn,gradIn.*spInfo.occEdge]);
mask= sum(TLM,3);
mask(find(~mask))=nan;
% figure(1); 
% subplot(1,2,1);imshow(-info.dis0,[]); title('before')
% subplot(1,2,2);imshow(-dpOut1,[]);title('after')
% [xx,yy] = meshgrid(1:size(dpOut1,1));
% figure(2);surf(xx,yy,dpOut1.*mask);

end

