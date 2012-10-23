function [nout,flagout]=nextn(u,v,n,flag)
if (u(n)<0.8*v(n)),
    if n<15 && flag==0,
        n=n+1;
    end
    flagout=flag;
    nout=n;
    return
end 
if n>7 && n<15 && flag==0,       
    n=n+1;
elseif n==15 &&flag==0,
    n=n-1;
    flag=1;
elseif flag==1 && n>7,
    n=n-1;
elseif flag==1 && n==7,
    flag=0;
end
flagout=flag;
nout=n;
end

