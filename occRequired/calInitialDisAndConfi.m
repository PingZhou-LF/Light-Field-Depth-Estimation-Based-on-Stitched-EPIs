function [InitialDis_Angle,confi]= calInitialDisAndConfi(cost)

[mRow,mCol,~]= size(cost);

[~,InitialDis_Angle]= min(cost,[],3); % figure;imshow(dpCorIdx,[]);
% [~,dpDefIdx]= min(defocus_response,[],3); % figure;imshow(dpDefIdx/101);

tmpVal= reshape(cost,mRow*mCol,[]);
[confVal1,~] = min(tmpVal,[],2);
minResp= reshape(confVal1,mRow,mCol); 
confVal2= mean(tmpVal,2);
meanResp= reshape(confVal2,mRow,mCol); 
corWeight= meanResp./minResp; % figure;imshow(corWeight/10);
nanIdx= find(isnan(corWeight));
corWeight(nanIdx)= 0.00001;

% tmpVal= reshape(defocus_response,mRow*mCol,[]);
% [confVal1,~] = min(tmpVal,[],2);
% confVal2= mean(tmpVal,2);
% defWeight= reshape(confVal2./confVal1,mRow,mCol); % figure;imshow(defWeight/10);

confMaxVal= max(min(max(corWeight(:)),5*median(corWeight(:)))); % ,...
% min(max(defWeight(:)),5*median(defWeight(:))));

confi= corWeight/confMaxVal;
confi(confi>1)= 1; % figure;imshow(cw);
% dw= defWeight/confMaxVal;
% dw(dw>1)=1; % figure;imshow(dw);
%%%%%%%加入可靠度随深度增加而降低的约束
% d_res = info.angle_res;
% du = ceil(d_res/2);
% mask = InitialDis_Angle >= du;
% maskdepth= mask.*InitialDis_Angle;
% mind = du;
% maxd = max(maskdepth(:));
% regmaskdepth = (maskdepth-mind)/(maxd-mind);
% regmaskdepth = mask.*regmaskdepth;
% weight = exp(-regmaskdepth);
% cw=cw.*weight;
%%%%%%%
end