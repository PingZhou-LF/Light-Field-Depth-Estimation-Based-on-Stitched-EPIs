function [cwn]= cwManipulate(spInfo)
spInfo.cwParams.lcwEdgeEnforce= 2;
spInfo.cwParams.Th_textureSupress= 0;

cw= spInfo.weight;

% 1. low confidence occlusion border shrinkage
dpDiff= spInfo.dis_initial-spInfo.dis_SP; % figure;imshow(dpDiff);
dpStretch= (dpDiff<0)*2./(1+exp(-dpDiff))+ (dpDiff>=0);
% figure;imshow(dpStretch);

% 2. high confidece intricate structure preservation
hcMask= cw>0.5; 
x= (~imopen(hcMask,strel('square',5))).*hcMask; % figure;imshow(x);
dpStretch= (~x).*dpStretch+ x; % figure;imshow(dpStretch);

% 3. high est. variation region supression
estVar = stdfilt(spInfo.dis_initial/spInfo.d_res,ones(3,3));
estVar= estVar/min(max(estVar(:)),10*median(estVar(:)));
Th_minEstVar= spInfo.cwParams.Th_minEstVar; % 0.1, 0.2
stimes= 10000; % stimes=10000;
estVarStretch= (estVar<=Th_minEstVar)+ (estVar>Th_minEstVar)*2./(1+exp(stimes*(estVar- Th_minEstVar))); 
% figure;imshow(estVarStretch);

cwn= cw.*dpStretch.*estVarStretch;
