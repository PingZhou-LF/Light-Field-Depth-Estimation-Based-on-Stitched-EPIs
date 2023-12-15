% 遮挡优化（1）——半拼接极图计算 

clear
addpath(genpath(pwd))
% nameAll={'boxes','antinous','greek','kitchen',...
%     'museum','pens','pillows','rosemary', ...
%     'table','tower','vinyl'};%

name='rosemary';
%%
pathLoad = 'D:\SLQ\codePackage\rosemary\data0\rosemary.mat';
load(pathLoad,'info');%info,dis0
pathLoad = 'D:\SLQ\codePackage\rosemary\dataOcc\rosemaryFinal.mat';
load(pathLoad);%info,dis0
pathLoad = 'D:\SLQ\codePackage\rosemary\normalization\rosemary.mat';
load(pathLoad,'LF');%
cfg.disp_max=3;
cfg.disp_min=-3;
%% 参342数
info.LF=LF;
[info.Nt,info.Ns,info.Ny,info.Nx,info.cc]=size(info.LF);
info.t0 = ceil(info.Nt/2);
info.s0 = ceil(info.Ns/2);
info.csai = squeeze(info.LF(info.t0,info.s0,:,:,:));
info.LF_Remap    = reshape(permute(info.LF, ...
    [1 3 2 4 5]), [info.Nt*info.Ny info.Ns*info.Nx 3]); % Remap图像
sai=cell(info.Nt,info.Ns);

dis_max=cfg.disp_max;
dis_min=cfg.disp_min;
info.dis_max=dis_max;
info.dis_min=dis_min;
info.alpha1=0.5;
info.alpha2=0.2;
info.res=181;
info.dpth=0.07;
info.angle_min = (atan(1/info.dis_max)*180)/pi; %弧度 转换成 °
info.angle_max = 180+(atan(1/info.dis_min)*180)/pi;
info.step = (info.angle_max-info.angle_min)/(info.res-1);
%%
for t_ind=1:1:info.Nt
    for s_ind=1:1:info.Ns
        sai{t_ind,s_ind} = squeeze(info.LF(t_ind,s_ind,:,:,:));
    end
end


%% SX
Ny=info.Ny;
Nx=info.Nx;
dis0=info.dis0;
flag_all=zeros(Ny,Nx);
for y_=1:1:Ny
    for x_=1:1:Nx
        xl= threshold(x_-4,Nx,1);
        xr= threshold(x_+4,Nx,1);
        dis_=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(x_-3,Nx,1):threshold(x_+3,Nx,1))));
        disl=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(xl-3,Nx,1):threshold(xl+3,Nx,1))));
        disr=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(xr-3,Nx,1):threshold(xr+3,Nx,1))));
        if dis_<disr && disl<disr %\
            flag_all(y_,x_)=1;
        else
            flag_all(y_,x_)=2;
        end

    end
end
% figure; imagesc(flag_all); axis equal; axis off;title('flag')

%%
% 
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
%%
tic
dis_step=0.1;
dis_half01_cpp=zeros(Ny,Nx);
dis_res=ceil((cfg.disp_max-cfg.disp_min)/dis_step + 1);
cost_half01_cpp = zeros(Ny, Nx, dis_res);
cost_half01_shape = zeros(Ny, Nx, dis_res);
%%
tic
[cost_half01_cpp, dis_half01_cpp] = disRefine_mex(sai, info.csai, flag_all, cost_half01_shape,dis_step);
toc
%% TY
temp=info.Nx;%%\
info.Nx=info.Ny;
info.Ny=temp;

sai=cell(info.Nt,info.Ns);
for t_ind=1:1:info.Nt
    for s_ind=1:1:info.Ns
        sai{s_ind,info.Nt+1-t_ind} = imrotate(squeeze(info.LF(t_ind,s_ind,:,:,:)),-90);
    end
end

% for i=1:1:5
%     figure(1);imshow(sai{5,i});pause(1)
% end

info.csai=imrotate(info.csai,-90);
dis0= imrotate(dis0,-90);

Ny=info.Ny;
Nx=info.Nx;
flag_all=zeros(Ny,Nx);
for y_=1:1:Ny
    for x_=1:1:Nx
        xl= threshold(x_-4,Nx,1);
        xr= threshold(x_+4,Nx,1);
        dis_=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(x_-3,Nx,1):threshold(x_+3,Nx,1))));
        disl=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(xl-3,Nx,1):threshold(xl+3,Nx,1))));
        disr=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(xr-3,Nx,1):threshold(xr+3,Nx,1))));
        if dis_<disr && disl<disr %涓?
            flag_all(y_,x_)=1;
        else
            flag_all(y_,x_)=2;
        end

    end
end
figure;imagesc(flag_all); axis equal; axis off;title('flag')



tic
dis_step=0.1;
dis_half02=zeros(Ny,Nx);
dis_res=ceil((info.dis_max-info.dis_min)/dis_step);
% cost_half02=zeros(info.Ny,info.Nx,dis_res);
for x=1:1:Nx
    for y=1:1:Ny
        [cost_half,dis_half02(y,x)] = disRefine(info.dis_max,info.dis_min,dis_step,info,x,y,flag_all(y,x),sai);
        cost_half02(y,x,:)=cost_half;
    end
end

toc

subplot(1,3,2);imshow(dis0,[],'InitialMagnification','fit');title('dis0')
subplot(1,3,3);imshow(dis_half02,[],'InitialMagnification','fit');title('dis_half')




path=['E:\0_research\1_code\lxy\',nameDataset,'\dataOcc\',name,'Occ'];
save(path,'dis_half01','cost_half01','dis_half02','cost_half02');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except ind nameAll nameDataset stSize

