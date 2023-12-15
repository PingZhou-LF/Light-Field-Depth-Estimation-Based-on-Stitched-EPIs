function cost = Initial_Depth_Compute(info)
%INITIAL_DEPTH_COMPUTE 此处显示有关此函数的摘要
%   此处显示详细说明
tic
[cost.slope_cost_sx,cost.gx_cost_sx,cost.gy_cost_sx]=mymex1(info.LF_Remap,info.Nx,info.Ny,info.Nt,...
                                      info.angle_min,info.angle_max,info.res);
%sx
% 转tyLF
toc
a='sx is done'

info.LF = permute(info.LF,[2 1 4 3 5]);
[info.Nt,info.Ns,info.Ny,info.Nx,info.c]=size(info.LF);
info.t0 = ceil(info.Nt/2);
info.s0 = ceil(info.Ns/2);
info.csai = squeeze(info.LF(info.t0,info.s0,:,:,:));
info.LF_Remap    = reshape(permute(info.LF, ...
[1 3 2 4 5]), [info.Nt*info.Ny info.Ns*info.Nx 3]); % Remap图像
tic
% 计算（xxx_cost是角度值，后缀sx表示用的是横向极图）
[cost.slope_cost_ty,cost.gx_cost_ty,cost.gy_cost_ty]=mymex1(info.LF_Remap,info.Nx,info.Ny,info.Nt,...
                                      info.angle_min,info.angle_max,info.res);
a='ty is done'
toc
cost.slope_cost_ty = permute(cost.slope_cost_ty,[2 1 3]);
cost.gx_cost_ty= permute(cost.gx_cost_ty,[2 1 3]);
cost.gy_cost_ty= permute(cost.gy_cost_ty,[2 1 3]);

end

