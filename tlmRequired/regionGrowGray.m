function [imgSelect,label] = regionGrowGray(flagPixel,x,y,img,threhold)

%  区域生长法分割彩色图像
f=img;  
[M,N]=size(f);
label=zeros(M,N); 
imgSelect=zeros(M,N);

y1=round(x);x1=round(y);
seed=f(x1,y1);       
 
label(x1,y1)=1;               %将Y中与所取点相对应位置的点设置为1
imgSelect(x1,y1)=f(x1,y1);
count=1;                  %记录每次判断一点周围八点符合条件的新点的数目

while count>0
count=0;
  for i=1:M
   for j=1:N
     if label(i,j)==1
       if (i-1)>0 && (i+1)<=M && (j-1)>0 && (j+1)<=N       %判断此点是否为图像边界上的点
         for u= -1:1                                        %判断点周围4点是否符合域值条件
           for v= -1:1                                       %u,v为偏移量             
              temp=f(i+u,j+v);
                %判断点是否未存在于生长区域中，并且为符合域值条件的点
             if  label(i+u,j+v)==0 && flagPixel(i+u,j+v)==0 && abs(temp-seed)<threhold
                 imgSelect(i+u,j+v)=f(i+u,j+v);
                 label(i+u,j+v)=1; 
                 flagPixel(i+u,j+v)=1;
                 count=1;
             end
           end  
         end
        if(count==0) 
            return;
        end
      end
     end
   end  
  end 
end




