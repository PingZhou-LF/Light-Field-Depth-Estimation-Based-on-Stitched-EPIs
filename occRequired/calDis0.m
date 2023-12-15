function [costHalf,dis_re] = calDis0(dis_max,dis_min,dis_step,info,x0,y0,sai)
dis_res=ceil((dis_max-dis_min)/dis_step);
cost=zeros(dis_res,2);

i=1;
%遍历角度计算cost
for dis_temp=dis_min:dis_step:dis_max
    slope=1./dis_temp;
    [~,~,cost_] = indexCalAll(info,slope,x0,y0,sai(:,:,:));
    cost(i,1)=cost_;
    cost(i,2)=slope;
    cost(i,3)=dis_temp;
    i=i+1;
end
b=cost(:,1);
[~,i]=min(b);
dis_re=cost(i,3);
costHalf=cost(:,1);

