function [index_x,index_y,cost_sum]= indexCalAll(info,k,x0,y0,sai)
% input：
%   (x0,y0):待优化点在深度图中的坐标
%   k:经过待优化点的某条直线的斜率（左下为原点
% output：
%   index：n*2,n=Nt*s0,某斜率k下待优化点在所有(s,t)内对应的(x,y),上->下
t0=info.t0;
s0=info.s0;
Nt=info.Nt;
Ns=info.Ns;
Nx=info.Nx;
Ny=info.Ny;
index_x=zeros(Ns*Nt,1);
index_y=zeros(Ns*Nt,1);
for delta_t=-(t0-1):1:(t0-1) %t循环，下-上
    index_x0  =  x0 + delta_t*Ns/k;%当前t对应的中心x
    index_y_   =  y0 + delta_t  /k;
    [index_x_temp,index_y_temp] = meshgrid(index_x0+(-(s0-1):1:(s0-1))/k,index_y_);
    index_x( (1+Ns*(delta_t+4) ):(Ns+Ns*(delta_t+4)))=index_x_temp;
    index_y( (1+Ns*(delta_t+4) ):(Ns+Ns*(delta_t+4)))=index_y_temp;
end

cost=zeros(Ns*Nt,1);
center=info.csai(y0,x0,:);%中心子孔径中当前点的color
for t_ind=1:1:Nt
    for s_ind=1:1:Ns
        x_ind  = index_x( (t_ind-1)*Ns+s_ind)-(t_ind-t0)*Ns/k;
        y_ind  = index_y( (t_ind-1)*Ns+s_ind);
        x_ceil = ceil (x_ind);x_ceil=threshold(x_ceil,Nx,1);
        x_floor= floor(x_ind);x_floor=threshold(x_floor,Nx,1);
        y_ceil = ceil (y_ind);y_ceil=threshold(y_ceil,Ny,1);
        y_floor= floor(y_ind);y_floor=threshold(y_floor,Ny,1);  
        dis_x_c = x_ceil-x_ind;
        dis_x_f = 1-dis_x_c;
        dis_y_c = y_ceil-y_ind;
        dis_y_f = 1-dis_y_c;     
        temp    = dis_y_c*dis_x_c*sai{-t_ind+Nt+1,-s_ind+Ns+1}(y_floor,x_floor,:) + ...
                  dis_y_c*dis_x_f*sai{-t_ind+Nt+1,-s_ind+Ns+1}(y_floor,x_ceil ,:) + ...
                  dis_y_f*dis_x_c*sai{-t_ind+Nt+1,-s_ind+Ns+1}(y_ceil ,x_floor,:) + ...
                  dis_y_f*dis_x_f*sai{-t_ind+Nt+1,-s_ind+Ns+1}(y_ceil ,x_ceil,: );
        cost_=sum(abs(temp-center))/3;
        cost((t_ind-1)*Ns+s_ind)=cost_*cost_;
    end
end

cost_sum=sum(cost);