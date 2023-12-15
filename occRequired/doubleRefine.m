function [doubleSX,disSX,doubleTY,disTY]=doubleRefine(dis_occ,info,flagOthers)
%% 
info.SPInnerNum = 15*15;
info.compactness=15;
[Rlabel,Esp]= SuperPixelSegment(info);
disSmooth =  medfilt2(dis_occ, [2 2],'symmetric');
disSmooth =  medfilt2(disSmooth, [2 2],'symmetric');
kSmooth = 1./disSmooth;
%%
tic
% doubleSX = doubleLinesRefine(info.LF_Remap,info.Nx,info.Ny,info.Nt,info.Ns...
%                                 ,kSmooth,double(flagOthers),81.0,double(Rlabel),1);
doubleSX = doubleLines1(info.LF_Remap,info.Nx,info.Ny,info.Nt,info.Ns...
                                 ,kSmooth,double(flagOthers),81.0,double(Rlabel),1);
toc

disSX=1./doubleSX;
disSX(disSX<info.dis_min)=info.dis_min;
disSX(disSX>info.dis_max)=info.dis_max;
%%
info.LF = permute(info.LF,[2 1 4 3 5]);
[info.Nt,info.Ns,info.Ny,info.Nx,info.c]=size(info.LF);
info.t0 = ceil(info.Nt/2);
info.s0 = ceil(info.Ns/2);
info.csai = squeeze(info.LF(info.t0,info.s0,:,:,:));
info.LF_Remap    = reshape(permute(info.LF, ...
[1 3 2 4 5]), [info.Nt*info.Ny info.Ns*info.Nx 3]); % Remap图像
tic
doubleTY = doubleLinesRefine(info.LF_Remap,info.Nx,info.Ny,info.Nt,info.Ns...
                                ,kSmooth',double(flagOthers'),81.0,double(Rlabel'),1); 
toc
disTY=1./doubleTY;
disTY(disTY<info.dis_min)=info.dis_min;
disTY(disTY>info.dis_max)=info.dis_max;
%%
% figure;
% subplot(1,2,1);imshow(disSX,[]);
% subplot(1,2,2);imshow(disTY',[]);
doubleTY=doubleTY';
disTY=disTY';