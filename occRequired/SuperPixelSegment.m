function [Rlabel,spBoundary] = SuperPixelSegment(info)
%SUPERPIXELSEGMENT 此处显示有关此函数的摘要
%   此处显示详细说明
numSP = info.SPInnerNum;
compactness = info.compactness;
csai = info.csai;
[Ny,Nx,~] = size(csai);
[Rlabel,~]= slicmex(im2uint8(csai),round(Ny*Nx/numSP),compactness);
    
% 标签化为1-K
K=length(unique(Rlabel)); %标签数
if K~=(range(Rlabel(:))+1)
    j=1;
  for i= min(Rlabel(:)):max(Rlabel(:))
        ind = find(Rlabel == i);
        if ~isempty(ind)
            Rlabel(ind)=j;
            j=j+1;
        end
  end
end

% 超像素分割边界
[Gx,Gy]= gradient(double(Rlabel));
spBoundary= ((Gx.^2+Gy.^2)~=0); 
end

