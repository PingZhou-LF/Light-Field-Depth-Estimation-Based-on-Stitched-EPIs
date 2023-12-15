function [costHalf,dis_re] = disRefine(dis_max,dis_min,dis_step,info,x0,y0,flag,sai)
dis_res=ceil((dis_max-dis_min)/dis_step);
cost=zeros(dis_res,2);

if flag==1%上 （后景在左前景在右
    i=1;
    %遍历角度计算cost
    for dis_temp=dis_min:dis_step:dis_max
       slope=1./dis_temp;
       [~,~,cost_] = indexCalUp(info,slope,x0,y0,sai(:,1:info.s0,:));
       cost(i,1)=cost_;
       cost(i,2)=slope;
       cost(i,3)=dis_temp;
       i=i+1;
    end
    b=cost(:,1);
    [~,i]=min(b);
    dis_re=cost(i,3);
    costHalf=cost(:,1);
    return
end


if flag==2%下 （后景在右前景在左
    i=1;
    %遍历角度计算cost
    for dis_temp=dis_min:dis_step:dis_max
       slope=1./dis_temp;
       [~,~,cost_] = indexCalDown(info,slope,x0,y0,sai(:,info.s0:info.Ns,:));
       cost(i,1)=cost_;
       cost(i,2)=slope;
       cost(i,3)=dis_temp;
       i=i+1;
    end
    b=cost(:,1);
    [~,i]=min(b);
    dis_re=cost(i,3);
    costHalf=cost(:,1);
    return
end