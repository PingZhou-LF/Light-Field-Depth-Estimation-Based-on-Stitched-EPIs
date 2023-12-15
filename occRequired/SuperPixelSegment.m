function [Rlabel,spBoundary] = SuperPixelSegment(info)
%SUPERPIXELSEGMENT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
numSP = info.SPInnerNum;
compactness = info.compactness;
csai = info.csai;
[Ny,Nx,~] = size(csai);
[Rlabel,~]= slicmex(im2uint8(csai),round(Ny*Nx/numSP),compactness);
    
% ��ǩ��Ϊ1-K
K=length(unique(Rlabel)); %��ǩ��
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

% �����طָ�߽�
[Gx,Gy]= gradient(double(Rlabel));
spBoundary= ((Gx.^2+Gy.^2)~=0); 
end

