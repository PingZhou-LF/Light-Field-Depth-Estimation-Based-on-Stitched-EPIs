function dis_filtered = disFilter(x_center,y_center,dis_re,info,w)

sigma = 0.01; 
csai=info.csai;
Nx=info.Nx;
Ny=info.Ny;
weight=0;%各点权重
weight_sum=0;
dis_sum=0;
var=0;
for i=-w:1:w
    for j=-w:1:w      
        x_temp=threshold(x_center + i,Nx,1);
        y_temp=threshold(y_center + j,Ny,1);
        var=sum(csai(y_temp,x_temp,:)-csai(y_center,x_center,:));
        weight=exp(-var*var/2*sigma*sigma);
        weight_sum=weight_sum+weight;
        dis_sum=dis_sum+weight*dis_re(y_temp,x_temp);
    end 
end


dis_filtered=dis_sum/weight_sum;
