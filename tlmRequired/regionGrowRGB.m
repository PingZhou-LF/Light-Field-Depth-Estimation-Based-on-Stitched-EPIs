function [imgSelect,label] = regionGrowRGB(flagPixel,x,y,img,threhold)

%  �����������ָ��ɫͼ��
f=img;   %monkey.jpg   clumsy2.bmp  za.bmp

[M,N]=size(f(:,:,1));
label=zeros(M,N); 
imgSelect=zeros(M,N,3);

y1=round(x);x1=round(y);
seed=[f(x1,y1,1),f(x1,y1,2),f(x1,y1,3)];       
 
label(x1,y1)=1;               %��Y������ȡ�����Ӧλ�õĵ�����Ϊ1
imgSelect(x1,y1,:)=f(x1,y1,:);
count=1;                  %��¼ÿ���ж�һ����Χ�˵�����������µ����Ŀ

while count>0
count=0;
  for i=1:M
   for j=1:N
     if label(i,j)==1
       if (i-1)>0 && (i+1)<=M && (j-1)>0 && (j+1)<=N       %�жϴ˵��Ƿ�Ϊͼ��߽��ϵĵ�
         for u= -1:1                                        %�жϵ���Χ4���Ƿ������ֵ����
           for v= -1:1                                       %u,vΪƫ����             
              r=f(i+u,j+v,1);g=f(i+u,j+v,2);b=f(i+u,j+v,3);
                %�жϵ��Ƿ�δ���������������У�����Ϊ������ֵ�����ĵ�
             if  label(i+u,j+v)==0 & abs(r-seed(1))<threhold && abs(g-seed(2))<threhold && abs(b-seed(3))<threhold && flagPixel(i+u,j+v)==0
                 imgSelect(i+u,j+v,:)=f(i+u,j+v,:);
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




