function [mout,nout,flag]=growth2d(m,n,u,v,sourcem,sourcen,t1,t2,flagin)
testm=sourcem+3;
testn=sourcen+3;
mout=m;
nout=n;
flag=flagin;
if flagin==0,
    if u((testm-1)*n+testn)<t1 && v((testm-1)*n+testn)>t2,
        flag=1;
    end
end
if flag==0,
    return
elseif flag==1 || flag==1.5,
    if ((m-sourcem)<(n-sourcen) && flag==1) || (m-sourcem+80<n-sourcen && flag==1.5),
        m=m+1;
        if m>80,
            m=1;
            flag=1.5;
        end
    elseif n==sourcen+40,
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
    elseif m==5+sourcem,
        flag=0;
    else
        n=n-1;
    end
elseif flag==2.5,
    if m-sourcem+80>n-sourcen,
        m=m-1;
        if m<1,
            m=80;
            flag=2;
        end
    else
        n=n-1;
    end
end
mout=m;
nout=n;
return
        
        


        