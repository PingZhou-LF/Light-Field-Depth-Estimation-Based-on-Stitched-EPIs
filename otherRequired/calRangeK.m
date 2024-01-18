function [kmin,kmax] = calRangeK(k0,b0,n)
if(k0>=1)%关于y对称
    k=1/k0;
    b=-b0/k0;   
    [n,q,p,s,k1,k2] = calRange(k,b,n);%1小2大
    kmin=1/k2;
    kmax=1/k1;
else
    if(k0>0)
        [~,~,~,~,k1,k2] = calRange(k0,b0,n);
        kmin=k1;
        kmax=k2;
    else
        if(k0>=-1)
            [~,~,~,~,k1,k2] = calRange(-k0,b0,n);
            kmin=-k2;
            kmax=-k1;
        else
                k=-1/k0;
                b=b0/k0;
                [~,~,~,~,k1,k2] = calRange(k,b,n);
                kmin=-1/k1;
                kmax=-1/k2;
         end
     end
end

if(isinf(kmax)) 
    kmax=k0;
end
if(isinf(kmin)) 
    kmin=k0;
end
