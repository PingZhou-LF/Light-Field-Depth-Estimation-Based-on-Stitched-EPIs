包含demoOcc中计算深度和代价的部分
原代码：
% tic
% dis_step=0.1;
% dis_res=ceil((info.dis_max-info.dis_min)/dis_step + 1);
% dis_half01=zeros(Ny,Nx);
% % cost_half01 = zeros(Ny, Nx, dis_res);
% 
% for x=1:1:Nx
%     for y=1:1:Ny
%         [cost_half,dis_half01(y,x)] = disRefine(info.dis_max,info.dis_min,dis_step,info,x,y,flag_all(y,x),sai);
%         cost_half01(y,x,:)=cost_half;
%     end
% end
% toc

改写后：
tic
[cost_half01_cpp, dis_half01_cpp] = disRefine_mex(sai, info.csai, flag_all, cost_half01_shape,dis_step);
toc