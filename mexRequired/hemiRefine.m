function [disLeftRefine,disRightRefine,disDownRefine,disUpRefine]=hemiRefine(LF,dis_occ,info,flagOcc,kRes)

[Nt,Ns,Ny,Nx,~]=size(LF);
maskLeft   =zeros(Ny,Nx);maskLeft(flagOcc==1) = 1;
maskRight = zeros(Ny,Nx);maskRight(flagOcc==2)= 1;
maskUp    = zeros(Ny,Nx);maskUp(flagOcc==4)   = 1;
maskDown  = zeros(Ny,Nx);maskDown(flagOcc==3) = 1;

sai=cell(Nt,Ns);
for t_ind=1:1:info.Nt
    for s_ind=1:1:info.Ns
        sai{t_ind,s_ind} = squeeze(info.LF(t_ind,s_ind,:,:,:));
    end
end
%% 右遮挡优化（后景在左前景在右
% tic
% disLeftRefine=LeftRefine(LF,dis_occ,info,maskLeft,kRes,sai);
% toc
%%
tic
disLeftRefine = LeftRefine_mex(Nt, Ns, dis_occ, info, maskLeft, kRes, sai);
toc
%% 左遮挡优化（后景在右前景在左
% disRightRefine=RightRefine(LF,dis_occ,info,maskRight,kRes,sai);
%%
tic
disRightRefine = RightRefine_mex(Nt, Ns, dis_occ, info, maskRight, kRes, sai);
toc
%%
temp=info.Nx;%%注意此处是和sx极图相反的
info.Nx=info.Ny;
info.Ny=temp;

sai=cell(info.Nt,info.Ns);
for t_ind=1:1:info.Nt
    for s_ind=1:1:info.Ns
        sai{s_ind,info.Nt+1-t_ind} = imrotate(squeeze(info.LF(t_ind,s_ind,:,:,:)),-90);
    end
end
info.csai=imrotate(info.csai,-90);
dis_occ= imrotate(dis_occ,-90);
LF= permute(LF,[2 1 4 3 5]);
%% 下遮挡优化
maskDown=imrotate(maskDown,-90);
disDownRefine=LeftRefine(LF,dis_occ,info,maskDown,kRes,sai);
disDownRefine=imrotate(disDownRefine,90);
%% 上遮挡优化
maskUp=imrotate(maskUp,-90);
disUpRefine=RightRefine(LF,dis_occ,info,maskUp,kRes,sai);
% figure;imshow(abs(dis_occ-disUpRefine),[])
disUpRefine=imrotate(disUpRefine,90);





