clearvars -except source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths

load(paths.dis0_path,'info');
load(paths.disOccCost_path);
load(paths.source_path,'LF');
dis0=info.dis0;
[Ny,Nx]=size(dis0);
csai=squeeze(LF(5,5,:,:,:));
info.csai=csai;
info.Ny=Ny;
info.Nx=Nx;
flag_all=zeros(Ny,Nx);
alpha=1;%
%

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

%%

thres1=0.1;
varDis=abs(dis_half-dis0);
[occ_y,occ_x]=find(varDis>thres1 );

w1=7;
hs1=5;
hr1=20;
deno1=8;
key_all=zeros(Ny,Nx);
w2=10;
hs2=5;
hr2=20;
deno2=8;
key_half=zeros(Ny,Nx);
dis_min=min(min(dis0));
dis_max=max(max(dis0));

for i=1:size(occ_x,1)
    x1= threshold(occ_x(i),Nx-w1,1+w1);
    y1= threshold(occ_y(i),Ny-w1,1+w1);
    x2= threshold(occ_x(i),Nx-w2,1+w2);
    y2= threshold(occ_y(i),Ny-w2,1+w2);
    key_all(y1,x1) =occRefine(y1,x1,w1,dis0      ,dis_min ,dis_max ,hs1,hr1,deno1);
    key_half(y2,x2)=occRefine(y2,x2,w2,dis_half  ,dis_min ,dis_max ,hs2,hr2,deno2);
end
% key_all = occRefine_mex1(w1, dis0, info, hs1, hr1, deno1, occ_x, occ_y);
% key_half = occRefine_mex1(w2, dis_half, info, hs2, hr2, deno2, occ_x, occ_y);
%%
[yDis,xDis]=find(key_half==1);
w3=3;
hs3=5;
hr3=15;
deno3=8;


csai_min=rgb2gray(min(min(info.csai)));
csai_max=rgb2gray(max(max(info.csai)));
key_csai=zeros(info.Ny,info.Nx);

for i=1:size(xDis,1)
    x3= threshold(xDis(i),info.Nx-w3,1+w3);
    y3= threshold(yDis(i),info.Ny-w3,1+w3);
    key_csai(y3,x3)=occRefine(y3,x3,w3,rgb2gray(info.csai),csai_min ,csai_max ,hs3,hr3,deno3);
end
% key_csai = occRefine_mex1(w3, rgb2gray(info.csai), info, hs3, hr3, deno3, xDis, yDis);
key_final=zeros(info.Ny,info.Nx);
key_final(key_all==1 | key_csai==1 )=1;
[occ_y_final,occ_x_final]=find(key_final==1);

%%
dis_occ=dis0;
w=1;
for i=1:size(occ_x_final,1)
    x=occ_x_final(i);
    y=occ_y_final(i);
    dis_occ(y,x)=disFilter(x,y,dis_half,info,w);
    %     dis_occ(y,x)=dis_half01(y,x);

end

%%
occ.alpha=alpha;
occ.thres1=thres1;
occ.w1=w1;occ.hs1=hs1;occ.hr1=hr1;occ.deno1=deno1;
occ.w2=w2;occ.hs2=hs2;occ.hr2=hr2;occ.deno2=deno2;
occ.w3=w3;occ.hs3=hs3;occ.hr3=hr3;occ.deno3=deno3;
occ.w=w;


save(paths.disOccFinal_path, 'occ','flag_all','dis_half01','cost_half01','dis_half02','cost_half02', ...
    'key_all','key_half','key_csai','key_final','dis0','dis_occ');