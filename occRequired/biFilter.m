function out = biFilter( in , img,wsize,pow)

%COSTFILTER 此处显示有关此函数的摘要
%   根据中心子孔径图像，对cost进行双边滤波

[Ny,Nx,res]   = size(in);
csai          = img;
out           = zeros(Ny,Nx,res);
% params
sigma         = 0.03; 
% 
for i =1:res
    for y=1:Ny
        for x=1:Nx
            cost_sum_tmp=0;
            weight_sum_tmp = 0;    
            for v = -wsize:wsize
                yv = y+v;
                if (yv <1) yv =1; end
                if (yv >Ny) yv = Ny; end
                for u=-wsize:wsize
                    xu = x+u;
                    if (xu < 1) xu = 1;  end
                    if (xu >Nx) xu = Nx; end
                        
                    % 计算滤波器权重
                    color_dif = (csai(y,x,1)-csai(yv,xu,1))^2+...
                                (csai(y,x,2)-csai(yv,xu,2))^2+...
                                (csai(y,x,3)-csai(yv,xu,3))^2;
                    weight    =  exp(-power(color_dif,pow)/(2*sigma*sigma));   
%                      weight    =  exp(-power(color_dif,pow));  
                    cost_sum_tmp = cost_sum_tmp + weight * in(yv,xu,i);
                    weight_sum_tmp = weight_sum_tmp + weight;
                end
            end  
            if (cost_sum_tmp~=0)
                out(y,x,i) = cost_sum_tmp/weight_sum_tmp;
            else
                out(y,x,i) = 1;
            end
        end
    end
end
end







