function [n,q,p,s,k1,k2] = calElements(Chain)

n=size(Chain,1);

for q=1:1:n
    flag=1;
    for i=1:1:n-q
        if(Chain(i+q)~=Chain(i))
            flag=0;
            break;
        end
    end  
    if(flag)%
        break;
    end
end

q=min(q,n);

p=0;
for i=1:1:q
    p=p+Chain(i);
end

for s=0:1:q-1
    flag=1;
    for i=1:1:q
        temp = floor(p * (i - s) / q) - floor(p * (i - s - 1) / q);
        if(Chain(i)~=temp)
           flag=0;
        end
    end
    if(flag)
        break;
    end
end



%l
for l = 0:1: q - 1
    if (mod(l * p + 1, q) == 0) 
        break;
    end
end

Fs = s;
Ls = s + floor((n - s) / q) * q;

Fsl = s + l - floor((s + l) / q) * q;
Lsl = s + l + floor((n - s - l) / q) * q;

k1 = p / q - 1.0 / (q * (Ls - Fsl));
k2 = p / q + 1.0 / (q * (Lsl - Fs));

end

