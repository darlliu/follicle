function fol1d7
%like fol1d6 but with Runge-Kutta time stepping
%for some reason blows up :(



clc
close all
global N M dx;
N=100;
% discritize N pts along x
x=linspace (0,1,N);
dx=1/N;
% distance vector
% solve in [0,1]



M=3000;
% time points

dt=0.005; 
%discritized time;


r1=[0; 0];
r2=[0; 0];

R1=zeros(2,M);
R2=R1;
%receptor bound matrix
%initial at 0
global Src Rec Prc Ends; 

Src=[10 11];
%the linear source


Rec=[40 41];
%receptors source

Prc=[];
%precortex

Ends=[50 51];
%the initial ends;

%note: you can change the length of these vectors but remember to change R
%accordingly

u=zeros(N,1);
%prepare temp conc vectors
u(Ends)=1E-5;
%initial conditions, change to whatever.
v=u;
w=v;
y=w;
U=zeros(N,M);
U(:,1)=u; 
V=U;
W=V;
Y=W;
%total concentration-time matrix of ligands


t=0:dt:dt*(M-1);
%time vector

%parameters:
global d D Rtot1 Rtot2 kon1 koff1 kon2 koff2 kdeg1 kdeg2 kg1 kg2...
    kgen1 kgen2 kr1 kr2 counter;

Rtot1=1E-3;
Rtot2=Rtot1;
%max receptor amount
kon1=5;
koff1=1;
kon2=5;
koff2=1;
kdeg1=1E-3;
kdeg2=1E-3;
kg1=5E-7;
kg2=5E-7;
kgen1=2E-5;
kgen2=2E-5;
kr1=140;
kr2=140;
counter=100;

d=1E-3;
%diffusion rate


D=MakeLaplacian1D(N);
% diffusion matrix

%rates
% a1=2;
% a2=2;
% b1=1;
% b2=1;
%coupling constants

H=zeros(M,1);

for i= 2: M,
    if counter==1,
        Ends=Ends+grow(r1,r2,kg1,kg2);
        if Ends(2)>62,
            Prc=Ends-2;
        else
            Prc=[];
        end
        if Ends(2)>=100,
            Ends=[99,100];
        end
        if Ends(1)<=50,
            Ends=[50,51];
        end
        counter=100;
    end
    %growth phase
    
    [u2,v2,w2,y2,r11,r22]=rkd4(u,v,w,y,r1,r2,0,0,0,0,0,0,0,dt);
    %k1
    [u3,v3,w3,y3,r111,r222]=rkd4(u,v,w,y,...
        r1,r2,u2/2,v2/2,w2/2,y2/2,r11/2,r22/2,dt/2,dt);
    %k2
    
    [u4,v4,w4,y4,r1111,r2222]=rkd4(u,v,w,y ...
    ,r1,r2,u3/2,v3/2,w3/2,y3/2,r111,r222,dt/2,dt);
    %k3
    
    [u5,v5,w5,y5,r11111,r22222]=rkd4(u,v,w,y, ...
        r1,r2,u4,v4,w4,y4,r1111,r2222,dt,dt);
   %k4
    
    U(:,i)=u+(1/6)*(u2+2*u3+2*u4+u5);
    V(:,i)=v+(1/6)*(v2+2*v3+2*v4+v5);
    W(:,i)=w+(1/6)*(w2+2*w3+2*w4+w5);
    Y(:,i)=y+(1/6)*(y2+2*y3+2*y4+y5);
    R1(:,i)=r1+(1/6)*(r11+2*r111+2*r1111+r11111);
    R2(:,i)=r2+(1/6)*(r22+2*r222+2*r2222+r22222);
    H(i)=Ends(1);
    u=U(:,i);
    v=V(:,i);
    w=W(:,i);
    y=Y(:,i);
    r1=R1(:,i);
    r2=R2(:,i);
end



%plotting routines
figure
subplot(1,2,1)
plot(t,sum(R1));
title('receptor bound overtime')
legend('BMP_L_R')
subplot(1,2,2)
plot(t,sum(R2));
title('receptor bound overtime')
legend('Wnt_L_R')
figure
plot(t,H);
title('growth vs time')
fig1=figure;
plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
legend('BMP','wnt','Noggin','Dkk')
windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
Movie=moviein(100,fig1,windowsize);
Movie(:,1)=getframe(fig1,windowsize);
frame=2;

for i=101:M/100:M,

    plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
    legend('BMP','wnt','Noggin','Dkk')
    %plot(x,U(:,i))
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;

end

% size(Results)

% size(U)



% size(Stem)
 movie(fig1, Movie, 100,8,windowsize);
 
return
% function [u3,v3,w3,y3,r111,r222]= ...
%     rkd4add(u,v,w,y,r1,r2,u3,v3,w3,y3,r111,r222,dt)
% 
% u3=(u+u3)*dt;
% v3=(v+v3)*dt;
% w3=(w+w3)*dt;
% y3=(y+y3)*dt;
% r111=(r1+r111)*dt;
% r222=(r2+r222)*dt;
% return
function [u2,v2,w2,y2,r11,r22]=rkd4(u,v,w,y,r1,r2,uu,vv,ww,yy,rr1,rr2,dt,dt2)
global d D Rtot1 Rtot2 kon1 koff1 kon2 koff2 kdeg1 kdeg2...
    kgen1 kgen2 kr1 kr2;
global Src Rec Prc Ends; 
global N dx;
    gu=rxn(u,r1,kon1,koff1,Rtot1,Rec,dt,N);
    gv=rxn(v,r2,kon2,koff2,Rtot2,Rec,dt,N);


    gu2=gen(kgen1,Src,N,dt);
    gv2=gen(kgen2,Ends,N,dt);
    gw2=gen(kgen1*0.9,Ends,N,dt);
    
    gy2=gen(kgen2*0.2,Src,N,dt);
    if isempty(Prc),
    else
%        gu2=gu2+gen(kgen1,Prc,N,dt);
%        gy2=gy2+gen(kgen2*0.6,Prc,N,dt);
    end
    gv2=0;
    gw2=0;
    
    [gu3,gw3]=rxn2(u,w,kr1,dt);
    [gv3,gy3]=rxn2(v,y,kr2,dt);
    fu=-kdeg1*u*dt+gu+gu2+gu3;
    fv=-kdeg1*v*dt+gv+gv2+gv3;
    fw=-kdeg1*w*dt+gw2+gw3;
    fy=-kdeg1*y*dt+gy2+gy3;
    u2=step(d,D,u,uu,fu,dx,dt,N)*dt2;
    v2=step(d,D,v,vv,fv,dx,dt,N)*dt2;
    w2=step(d,D,w,ww,fw,dx,dt,N)*dt2;
    y2=step(d,D,y,yy,fy,dx,dt,N)*dt2;

    r11=step2(u,r1,rr1,Rtot1,kon1,koff1,kdeg1,Rec,dt)*dt2;
    r22=step2(v,r2,rr2,Rtot2,kon2,koff2,kdeg2,Rec,dt)*dt2;


return

function n=grow(r1,r2,kg11,kg22)
%this part is tricky. must return n integer and n<N
if sum(r1)<kg1 && sum(r2)>kg2,
    n=1;
elseif sum(r1)>kg11*1.5 && sum(r2)<kg22,
    n=-1;
else
    n=0;
end

return
function gr=rxn(u,r,kon,koff,rtot,Rec,dt,N)
%gives reaction part given previous quantities and rates
gr=zeros(N,1);
gr(Rec)=dt*(-kon*u(Rec).*(rtot-r)+koff*(r));

return

function gu2=gen(kgen1,Src,N,dt)
gu2=zeros(N,1);
gu2(Src)=dt*kgen1;

%gu2(Src)=dt*(sum(r1)/(kgen1+sum(r1)));
%gv2(Ends)=dt*(sum(r2)/(kgen2+sum(r2)));
return


function u2=step(d,D,u,uu,F,dx,dt,N)
%I=eye(N,N);
%R=zeros(N,1);
u2=zeros(size(u));
u2(1:N-1) =d*(dt/dx^2)*D*u(1:N-1)+F(1:N-1);
u2=u+uu*dt+u2;
%u =(I-d*(dt/dx^2)*D)\(u+F);
return

function [gv3,gy3]=rxn2(v,y,kr2,dt)
gv3=-v.*y*kr2*dt;
gy3=gv3;


function r=step2(u,r,rr,rtot,kon,koff,kdeg,Rec,dt)
r=(kon*u(Rec).*(rtot-r)-(koff+kdeg)*r)*dt+rr*dt+r;
return
% function ru=step3(u,r,rtot,kon,koff,Rec,dt)
% ru=u*0;
% ru(Rec)=-(kon*u(Rec).*(rtot-r)-(koff))*dt;
% return