clc;
clearvars -except source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths

% Argument Setting
load(paths.source_path,'LF','cfg');
info.LF=LF;
[info.Nt,info.Ns,info.Ny,info.Nx,info.cc]=size(info.LF);
info.t0 = ceil(info.Nt/2);
info.s0 = ceil(info.Ns/2);
info.csai = squeeze(info.LF(info.t0,info.s0,:,:,:));
info.LF_Remap    = reshape(permute(info.LF, ...
    [1 3 2 4 5]), [info.Nt*info.Ny info.Ns*info.Nx 3]); 

dis_max=cfg.disp_max;
dis_min=cfg.disp_min;
info.dis_max=dis_max;
info.dis_min=dis_min;
info.alpha1=0.5;
info.alpha2=0.2;
info.res=181;
info.dpth=0.07;
info.angle_min = (atan(1/info.dis_max)*180)/pi; 
info.angle_max = 180+(atan(1/info.dis_min)*180)/pi;
info.step = (info.angle_max-info.angle_min)/(info.res-1);

%% Initial Depth
costinfo = Initial_Depth_Compute(info);
%%
alpha1=info.alpha1;
alpha2=info.alpha2;
slope_cost = alpha1*costinfo.slope_cost_sx + (1-alpha1)*costinfo.slope_cost_ty;
gx_cost = alpha1*costinfo.gx_cost_sx + (1-alpha1)*costinfo.gx_cost_ty;
gy_cost = alpha1*costinfo.gy_cost_sx + (1-alpha1)*costinfo.gy_cost_ty;

combined_cost = (1-alpha2)*slope_cost + alpha2*(2*gy_cost);
[angle_ind,c]= calInitialDisAndConfi(combined_cost);
c = double((c-min(c(:)))/range(c(:)));
dis0 = 1./tan((info.angle_min + (angle_ind-1)*info.step)*pi/180);
dis0(dis0<info.dis_min)=info.dis_min;
dis0(dis0>info.dis_max)=info.dis_max;
figure(1);imshow(dis0,[])
info.dis0=dis0;
info.c = c;

save(paths.dis0_path,'info');