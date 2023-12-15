clearvars -except source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths

load(paths.source_path);
load(paths.disOccFinal_path);
load(paths.dis0_path,'info');

kRes=81;
info.LF=LF;
[info.Nt,info.Ns,info.Ny,info.Nx,info.cc]=size(info.LF);
info.t0 = ceil(info.Nt/2);
info.s0 = ceil(info.Ns/2);
info.csai = squeeze(info.LF(info.t0,info.s0,:,:,:));
info.LF_Remap    = reshape(permute(info.LF, ...
    [1 3 2 4 5]), [info.Nt*info.Ny info.Ns*info.Nx 3]); 
info.dis_min = cfg.disp_min;
info.dis_max = cfg.disp_max;
%% mask
[Ny,Nx]=size(flag_all);
flagOthers = zeros(Ny,Nx);flagOthers(key_final~=1) = 1;
flagOcc = flag_all;flagOcc(key_final~=1) = 0;
maskLeft   =zeros(Ny,Nx);maskLeft(flagOcc==1) = 1;
maskRight = zeros(Ny,Nx);maskRight(flagOcc==2)= 1;
maskUp    = zeros(Ny,Nx);maskUp(flagOcc==4)   = 1;
maskDown  = zeros(Ny,Nx);maskDown(flagOcc==3) = 1;
%%
[doubleSX,disSX,doubleTY,disTY]=doubleRefine(dis_occ,info,flagOthers);
%%
[disLeftRefine,disRightRefine,disDownRefine,disUpRefine]=hemiRefine(LF,dis_occ,info,flagOcc,kRes);
%% 

disLinesRefine = (disSX+disTY)/2.*(flagOthers) + disLeftRefine.* (maskLeft) + ...
                                                 disRightRefine.*(maskRight)+ ...
                                                 disUpRefine.*   (maskUp)   + ...
                                                 disDownRefine.* (maskDown);

dis1=disLinesRefine;

alpha=0.5;
for i=1:1:size(cost_half02,3)
    cost_half02_(:,:,i)=imrotate(cost_half02(:,:,i),90);
end
cost_half02=cost_half02_;
combined_cost = (1-alpha)*cost_half01+alpha*(2*cost_half02);
[angle_ind,c_occ]= calInitialDisAndConfi(combined_cost);
c0=info.c;
c_occ=c_occ.*(~flagOthers)+c0.*(flagOthers);
c_occ = double((c_occ-min(c_occ(:)))/range(c_occ(:)));


save(paths.disMLR_path, 'dis0','dis_occ','dis1','c0','c_occ');