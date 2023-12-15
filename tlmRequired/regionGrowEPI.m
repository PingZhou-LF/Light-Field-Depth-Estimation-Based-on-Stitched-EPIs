function [h,tlmLine,label,xl,xr] = regionGrowEPI(x0,y0,LF,threhold,dis,thresDis,type)
% 输出：
% h：Ns*Nx*3 区域RGB图像
% tlmLine：1*N，当前区域的x坐标
% label：Ns*Nx,区域RGB图像不为零的label=1；其余label=0
% xl:当前线段的xmin
% xr:当前线段的xmax
[Nt,Ns,Ny,Nx,~]=size(LF);
if type=='sx'
f=squeeze(LF((Nt+1)/2,:,y0,:,:));
x=(Ns+1)/2;
y=x0;
[M,N]=size(f(:,:,1));
label=zeros(M,N); 
h=zeros(M,N,3);
x1=round(x);y1=round(y);
seed=[f(x1,y1,1),f(x1,y1,2),f(x1,y1,3),dis(y0,x0)];       
% seed=[f(x1,y1,1),f(x1,y1,2),f(x1,y1,3)];       

label(x1,y1)=1;               %将Y中与所取点相对应位置的点设置为1
h(x1,y1,:)=f(x1,y1,:);
count=1;                  %记录每次判断一点周围八点符合条件的新点的数目

while count>0
count=0;
  for i=1:M
   for j=1:N
     if label(i,j)==1
       if (i-1)>0 && (i+1)<=M && (j-1)>0 && (j+1)<=N && y0-1>0 && y0+1<Ny && x0-1>0 &&  x0+1<Nx     %判断此点是否为图像边界上的点
         for u= -1:1                                        %判断点周围4点是否符合域值条件
           for v= -1:1                                       %u,v为偏移量             
              r=f(i+u,j+v,1);g=f(i+u,j+v,2);b=f(i+u,j+v,3);disTemp=dis(y0,x0+v);
                %判断点是否未存在于生长区域中，并且为符合域值条件的点
             if  label(i+u,j+v)==0 && abs(r-seed(1))<threhold && abs(g-seed(2))<threhold && abs(b-seed(3))<threhold ...
                     && abs(disTemp-seed(4))<thresDis
                 h(i+u,j+v,:)=f(i+u,j+v,:);
                 label(i+u,j+v)=1; 
                 count=1;
             end
           end  
        end
      end
     end
   end  
  end 
end
% figure;imshow(label,[],'InitialMagnification','fit')
tlmLine=find(h((Ns+1)/2,:,1)~=0);
% min(tlmLine)
% max(tlmLine)
x1_l=find(h(1,:,1)~=0, 1 );
x2_l=find(h(Ns,:,1)~=0, 1 );
x1_r=find(h(1,:,1)~=0, 1, 'last' );
x2_r=find(h(Ns,:,1)~=0, 1, 'last' );
xl=ceil(((Ns+1)/2-1)/(Ns-1)*(x2_l-x1_l)+x1_l);
xr=ceil(((Ns+1)/2-1)/(Ns-1)*(x2_r-x1_r)+x1_r);

end

if type=='ty'
f=squeeze(LF(:,(Ns+1)/2,:,x0,:));
x=(Nt+1)/2;
y=y0;
[M,N]=size(f(:,:,1));
label=zeros(M,N); 
h=zeros(M,N,3);
x1=round(x);y1=round(y);
seed=[f(x1,y1,1),f(x1,y1,2),f(x1,y1,3),dis(y0,x0)];       
% seed=[f(x1,y1,1),f(x1,y1,2),f(x1,y1,3)];       

label(x1,y1)=1;               %将Y中与所取点相对应位置的点设置为1
h(x1,y1,:)=f(x1,y1,:);
count=1;                  %记录每次判断一点周围八点符合条件的新点的数目

while count>0
count=0;
  for i=1:M
   for j=1:N
     if label(i,j)==1
       if (i-1)>0 && (i+1)<=M && (j-1)>0 && (j+1)<=N && y0-1>0 && y0+1<Ny && x0-1>0 &&  x0+1<Nx      %判断此点是否为图像边界上的点
         for u= -1:1                                        %判断点周围4点是否符合域值条件
           for v= -1:1                                       %u,v为偏移量             
              r=f(i+u,j+v,1);g=f(i+u,j+v,2);b=f(i+u,j+v,3);disTemp=dis(y0+u,x0);
                %判断点是否未存在于生长区域中，并且为符合域值条件的点
             if  label(i+u,j+v)==0 && abs(r-seed(1))<threhold && abs(g-seed(2))<threhold && abs(b-seed(3))<threhold  ...
                     && abs(disTemp-seed(4))<thresDis
                 h(i+u,j+v,:)=f(i+u,j+v,:);
                 label(i+u,j+v)=1; 
                 count=1;
             end
           end  
        end
      end
     end
   end  
  end 
end
% figure;imshow(label,[],'InitialMagnification','fit')
tlmLine=find(h((Ns+1)/2,:,1)~=0);
% min(tlmLine)
% max(tlmLine)
x1_l=find(h(1,:,1)~=0, 1 );
x2_l=find(h(Ns,:,1)~=0, 1 );
x1_r=find(h(1,:,1)~=0, 1, 'last' );
x2_r=find(h(Ns,:,1)~=0, 1, 'last' );
xl=ceil(((Ns+1)/2-1)/(Ns-1)*(x2_l-x1_l)+x1_l);
xr=ceil(((Ns+1)/2-1)/(Ns-1)*(x2_r-x1_r)+x1_r);
end
