function [y,x]=centerSelect(LF,areaAll)
[Nt,Ns,~,~,~]=size(LF);
csai=squeeze(LF((Nt+1)/2,(Ns+1)/2,:,:,:));
[height,width]=size(areaAll);
[L,num]=bwlabel(areaAll,4);
plot_x=zeros(1,num);%%用于记录质心位置的坐标
plot_y=zeros(1,num);

for k=1:num  %%num个区域依次统计质心位置
    sum_x=0;sum_y=0;area=0;
    for i=1:height
    for j=1:width
       if L(i,j)==k
        sum_x=sum_x+i;
        sum_y=sum_y+j;
        area=area+1;   
       end
    end
    end
    plot_x(k)=fix(sum_x/area);
    plot_y(k)=fix(sum_y/area);
end

figure;
imshow(imoverlay(csai,areaAll,'y'),[],'InitialMagnification','fit');
for i=1:num
hold on
plot(plot_y(i) ,plot_x(i), 'r+','MarkerSize',5);
end

y=plot_y';
x=plot_x';