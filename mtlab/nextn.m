function [nout,flagout,cnt]=nextn(u,v,n,flag,counter)
cnt=counter-1;
if cnt>1,
    nout=n;
    flagout=flag;
    return
else
    cnt=1000;
end
if flag==0 && n==7 && (u(7)<1.2E-3 && v(7)>1.9E-3),
    n=n+1;
    flagout=flag;
    nout=n;
    return
end 
if n>7 && n<15 && flag==0,       
    n=n+1;
elseif n==15 && flag==0,
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

