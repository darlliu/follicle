function [mout,nout,flag]=growth2d2(M,N,m,n,r1,r2,sourcem,sourcen,t1,t2,flagin)
mout=m;
nout=n;
flag=flagin;
if flagin==0,
    if r1<t1 && r2>t2,
        flag=1;
    end
end
if flag==0,
    return
elseif flag==1,
    if (m-sourcem)<(n-sourcen),
        m=m+1;
        if m>M,
            m=1;
            flag=1.5;
        end
    elseif n==sourcen+M/2,
        n=n-1;
        flag=flag+1;
        nout=n;
        return;
    else
        n=n+1;
    end
elseif flag==1.5,
    if (m-sourcem+M<n-sourcen ),
        m=m+1;
    elseif n==sourcen+M/2,
        n=n-1;
        flag=flag+1;
        nout=n;
        return;
    else
        n=n+1;
    end
elseif flag==2,
    if m-sourcem>n-sourcen,
        m=m-1;
    elseif m==3+sourcem,
        flag=0;
    else
        n=n-1;
    end
elseif flag==2.5,
    if m-sourcem+M>n-sourcen,
        m=m-1;
        if m<1,
            m=M;
            flag=2;
        end
    else
        n=n-1;
    end
end
mout=m;
nout=n;
return