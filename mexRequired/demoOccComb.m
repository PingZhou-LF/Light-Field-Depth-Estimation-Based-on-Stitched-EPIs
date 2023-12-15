% 遮挡优化（2）——遮挡区域计算+优化
clear all
addpath(genpath(pwd));
% name='boardgames';
% pathLoad=['E:\0_research\1_code\lxy\add\dataOcc\',name,'Occ'];
% load(pathLoad);
% pathLoad=['E:\0_research\1_code\lxy\add\',name];
% load(pathLoad);
% Ny=info.Ny;
% Nx=info.Nx;
% dis0=info.dis0;
% csai=info.csai;

clear
addpath(genpath(pwd))
% nameAll={'boxes','greek','kitchen',...
%     'pens','pillows','rosemary', ...
%     'table','tower','vinyl'};%

name = 'rosemary';
pathLoad = 'D:\SLQ\codePackage\rosemary\data0\rosemary.mat';
load(pathLoad,'info');%info,dis0
pathLoad = 'D:\SLQ\codePackage\rosemary\dataOcc\rosemaryFinal.mat';
load(pathLoad);%info,dis0
pathLoad = 'D:\SLQ\codePackage\rosemary\normalization\rosemary.mat';
load(pathLoad,'LF');%
dis0=info.dis0;
[Ny,Nx]=size(dis0);
csai=squeeze(LF(5,5,:,:,:));
info.csai=csai;
info.Ny=Ny;
info.Nx=Nx;
flag_all=zeros(Ny,Nx);
alpha=1;%
% 
tic
for y_=1:1:Ny
    for x_=1:1:Nx
        xl= threshold(x_-4,Nx,1);%
        xr= threshold(x_+4,Nx,1);%
        yu= threshold(y_-4,Ny,1);%
        yd= threshold(y_+4,Ny,1);%
        dis_=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(x_-3,Nx,1):threshold(x_+3,Nx,1))));
        
        disl=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(xl-3,Nx,1):threshold(xl+3,Nx,1))));
        disr=sum(sum(dis0(threshold(y_-3,Ny,1):threshold(y_+3,Ny,1),threshold(xr-3,Nx,1):threshold(xr+3,Nx,1))));
        disu=sum(sum(dis0(threshold(yu-3,Ny,1):threshold(yu+3,Ny,1),threshold(x_-3,Nx,1):threshold(x_+3,Nx,1))));
        disd=sum(sum(dis0(threshold(yd-3,Ny,1):threshold(yd+3,Ny,1),threshold(x_-3,Nx,1):threshold(x_+3,Nx,1))));
        if dis_<disr && disl<disr %
            right = abs(disr-dis_)+abs(disr-disl);
            up    = abs(disu-dis_)+abs(disu-disd);
            down  = abs(disd-dis_)+abs(disd-disu);
            if(dis_<disu && disd<disu)%
                if(up>right*alpha) 
                    flag_all(y_,x_)=3;
                else
                    flag_all(y_,x_)=1;
                end
            else
                if(dis_<disd && disu<disd)%
                    if(down>right*alpha) 
                        flag_all(y_,x_)=4;
                    else
                        flag_all(y_,x_)=1;
                    end
                else
                   flag_all(y_,x_)=1; 
                end
            end
        else%
            left  = abs(disl-dis_)+abs(disl-disr);
            up    = abs(disu-dis_)+abs(disu-disd);
            down  = abs(disd-dis_)+abs(disd-disu);
            if(dis_<disu && disd<disu)%
                if(up>left*alpha)  
                    flag_all(y_,x_)=3;
                else
                    flag_all(y_,x_)=2;
                end
            else
                if(dis_<disd && disu<disd)%
                    if(down>left*alpha)  
                        flag_all(y_,x_)=4;
                    else
                        flag_all(y_,x_)=2;
                    end
                else
                   flag_all(y_,x_)=2; 
                end
            end
        end 
    end
end
% figure; imagesc(flag_all); axis equal; axis off;
% figure; imshow(labeloverlay(info.csai,flag_all)); axis equal; axis off;title('flag')
toc

%% 
dis_half02_=imrotate(dis_half02,90);
dis_half=zeros(Ny,Nx);
[flagLRy,flagLRx]=find(flag_all<3);
[flagUDy,flagUDx]=find(flag_all>2);
for i=1:1:size(flagLRx,1)
    dis_half(flagLRy(i),flagLRx(i))=dis_half01(flagLRy(i),flagLRx(i));
end

for i=1:1:size(flagUDx,1)
    dis_half(flagUDy(i),flagUDx(i))=dis_half02_(flagUDy(i),flagUDx(i));
end
% figure;imshow(dis_half,[],'InitialMagnification','fit');
% figure;imshow(dis_half01,[],'InitialMagnification','fit');title('sx')
% figure;imshow(imrotate(dis_half02,90),[],'InitialMagnification','fit');title('ty')
% 
%  [InitialDis_Angle,confi]= calInitialDisAndConfi(cost_half01);figure;imshow(confi)
%  [InitialDis_Angle,confi]= calInitialDisAndConfi(cost_half02);figure;imshow(imrotate(confi,90))
%% 

thres1=0.1;
varDis=abs(dis_half-dis0);
[occ_y,occ_x]=find(varDis>thres1 );
% figure(1);imshow(info.csai,[],'InitialMagnification','fit');hold on
% plot(occ_x,occ_y,'r.','MarkerSize',2);
% figure;imagesc(varDis); axis equal;axis off;

w1=3;
hs1=5;
hr1=10;
deno1=8;
key_all=zeros(info.Ny,info.Nx);
w2=3;
hs2=5;
hr2=10;
deno2=8;
key_half=zeros(info.Ny,info.Nx);
dis_min=min(min(dis0));
dis_max=max(max(dis0));

%%
tic
for i=1:size(occ_x,1)
    x1= threshold(occ_x(i),info.Nx-w1,1+w1);
    y1= threshold(occ_y(i),info.Ny-w1,1+w1);
    x2= threshold(occ_x(i),info.Nx-w2,1+w2);
    y2= threshold(occ_y(i),info.Ny-w2,1+w2);  
    key_all(y1,x1) =occRefine(y1,x1,w1,dis0      ,dis_min ,dis_max ,hs1,hr1,deno1);
    key_half(y2,x2)=occRefine(y2,x2,w2,dis_half  ,dis_min ,dis_max ,hs2,hr2,deno2);
end
toc
% figure;imagesc(key_all);axis equal; axis off;title("maskAll")
% figure(1);imshow(imoverlay(csai,key_all),[],'InitialMagnification','fit');

% figure;imagesc(key_half);axis equal; axis off;title("maskHalf")
% figure;imshow(imoverlay(csai,key_half),[],'InitialMagnification','fit');
% key_dis=zeros(info.Ny,info.Nx);
% key_dis(key_all==1 | key_half==1 )=1;
% figure;imagesc(key_dis);axis equal; axis off;title("key_dis")
%%
% key_all_cpp = zeros(info.Ny,info.Nx);
% key_half_cpp = zeros(info.Ny,info.Nx);
tic
key_all_cpp = occRefine_mex1(w1, dis0, info, hs1, hr1, deno1, occ_x, occ_y);
toc
%%
tic
key_half_cpp = occRefine_mex1(w2, dis_half, info, hs2, hr2, deno2, occ_x, occ_y);
toc
%% 
[yDis,xDis]=find(key_half==1);
w3=6;
hs3=5;
hr3=15;
deno3=8;


csai_min=rgb2gray(min(min(info.csai)));
csai_max=rgb2gray(max(max(info.csai)));
% key_csai=zeros(info.Ny,info.Nx);
% tic
% for i=1:size(xDis,1)
%     x3= threshold(xDis(i),info.Nx-w3,1+w3);
%     y3= threshold(yDis(i),info.Ny-w3,1+w3); 
%     key_csai(y3,x3)=occRefine(y3,x3,w3,rgb2gray(info.csai),csai_min ,csai_max ,hs3,hr3,deno3);
% end
% toc

key_csai = occRefine_mex1(w3, rgb2gray(info.csai), info, hs3, hr3, deno3, xDis, yDis);

% figure(1); imagesc(key_csai);axis equal; axis off;title("maskCsai")
key_final=zeros(info.Ny,info.Nx);
key_final(key_all==1 | key_csai==1 )=1;
[occ_y_final,occ_x_final]=find(key_final==1);
figure(1);imshow(imoverlay(csai,key_final),[],'InitialMagnification','fit');
% figure; imagesc(key_final);axis equal; axis off;title("maskFinal")
% figure;imshow(labeloverlay(dis0,3*key_final),[],'InitialMagnification','fit');title("maskFinal")
%% 
tic
dis_occ=dis0;
w=1;
for i=1:size(occ_x_final,1)
    x=occ_x_final(i);
    y=occ_y_final(i);
    dis_occ(y,x)=disFilter(x,y,dis_half,info,w);
%     dis_occ(y,x)=dis_half01(y,x);
    
end
toc
% figure;imshow(dis0,[],'InitialMagnification','fit');title('disAll')
% figure;imshow(dis_half,[],'InitialMagnification','fit');title('disHalf')
figure(1);imshow(dis_occ,[],'InitialMagnification','fit');title('disOcc')
%% 
% 
% errorDis0=abs(gt-dis0);figure;imshow(errorDis0,[],'InitialMagnification','fit');title('errorDis0')
% errorDisOcc=abs(gt-dis_occ);figure;imshow(errorDisOcc,[],'InitialMagnification','fit');title('errorDisOcc')
% 
% rmseDis0      = sqrt(mean(mean((errorDis0).^2)))
% rmseDisDisOcc = sqrt(mean(mean((errorDisOcc).^2)))
% 
% [occ_y_final,occ_x_final]=find(key_final==1);
% numOcc=size(occ_x_final,1);
% rmseOcc0 = sqrt(sum(errorDis0(key_final==1))/numOcc)
% rmseOcc  = sqrt(sum(errorDisOcc(key_final==1))/numOcc)

occ.alpha=alpha;
occ.thres1=thres1;
occ.w1=w1;occ.hs1=hs1;occ.hr1=hr1;occ.deno1=deno1;
occ.w2=w2;occ.hs2=hs2;occ.hr2=hr2;occ.deno2=deno2;
occ.w3=w3;occ.hs3=hs3;occ.hr3=hr3;occ.deno3=deno3;
occ.w=w;

pathResult=['E:\0_research\1_code\lxy\',nameDataset,'\dataOcc\',name,'Final'];
% save(pathResult, 'occ','flag_all','dis_half01','cost_half01','dis_half02','cost_half02', ...%%%%%%%%%%%%%%%%%%%%%%%%%%
%     'key_all','key_half','key_csai','key_final','dis0','dis_occ','gt');
save(pathResult, 'occ','flag_all','dis_half01','cost_half01','dis_half02','cost_half02', ...%%%%%%%%%%%%%%%%%%%%%%%%%%
    'key_all','key_half','key_csai','key_final','dis0','dis_occ');
