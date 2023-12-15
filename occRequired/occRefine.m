function [key]=occRefine(y,x,w,dis0,minPatch,maxPatch,hs,hr,deno)
patch=dis0(y-w:y+w,x-w:x+w);
% figure;imshow(v,[],'InitialMagnification','fit')
v=255./(maxPatch-minPatch).*patch+255*minPatch/(minPatch-maxPatch);
l = meanshsegm(v,hs,hr) ;
%  figure;imshow(label2rgb(l-1),'InitialMagnification','fit') ; 
sort_initial=tabulate(l(:));
sort_result=sortrows(sort_initial,2);
if size(sort_result,1)==1
    key=0;
     else if sort_result(end-1,2)>size(patch,1)*size(patch,2)/deno
        key=1;
         else
        key=0;
     end
end