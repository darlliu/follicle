function fol1d6
%another attempt at the older model.

clc
close all
N=200;
% discritize N pts along x
x=linspace (0,1,N);
dx=1/N;
% distance vector
% solve in [0,1]



M=10000;
% time points

dt=0.5; 
%discritized time;


r1=[0; 0];
r2=[0; 0];

R1=zeros(2,M);
R2=R1;
%receptor bound matrix
%initial at 0


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

D=MakeLaplacian1D(N);
% diffusion matrix


t=0:dt:dt*(M-1);
%time vector

%parameters:

Rtot1=1E-3;
Rtot2=Rtot1;
%max receptor amount
kon1=5E-2;
koff1=1E-2;
kon2=5E-2;
koff2=1E-2;
kdeg1=1E-3;
kdeg2=1E-3;
kg1=1E-7;
kg2=4E-7;
kgen1=2.1E-5;
kgen2=2E-5;
kr1=140;
kr2=140;
flag=0;
%rates
% a1=2;
% a2=2;
% b1=1;
% b2=1;
%coupling constants
d=1E-3;
%diffusion rate
H=zeros(M,1);
counter=100;
for i= 2: M,
    gu=rxn(u,r1,kon1,koff1,Rtot1,Rec,dt,N);
    gv=rxn(v,r2,kon2,koff2,Rtot2,Rec,dt,N);
    if counter==1,
        [n,flag]=grow(r1,r2,kg1,kg2,flag);
        Ends=Ends+n;
        if Ends(2)>75,
            Prc=Ends-2;
        else
            Prc=[];
        end
        if Ends(2)>=100,
            Ends=[99,100];
        end
        if Ends(1)<=50,
            Ends=[50,51];
            flag=0;
        end
        counter=100;
    end
    counter=counter-1;
    H(i)=Ends(1);
    gu2=gen(kgen1*0.2,Src,N,dt);
    gv2=gen(kgen2*0.3,Ends,N,dt);
    gv2=gv2+gen(kgen2*0.3,Src,N,dt);
    gw2=gen(kgen1*0.9,Ends,N,dt);
    gy2=gen(kgen2*0.2,Src,N,dt);
    if isempty(Prc),
    else
        gu2=gu2+gen(kgen1*0.5,Prc,N,dt);
        gy2=gy2+gen(kgen2*0.9,Prc,N,dt);
    end
    [gu3,gw3]=rxn2(u,w,kr1,dt);
    [gv3,gy3]=rxn2(v,y,kr2,dt);
    
    U(:,i)=step(d,D,u,-kdeg1*u*dt+gu+gu2+gu3,dx,dt,N);
    V(:,i)=step(d,D,v,-kdeg1*v*dt+gv+gv2+gv3,dx,dt,N);
    W(:,i)=step(d,D,w,-kdeg1*w*dt+gw2+gw3,dx,dt,N);
    Y(:,i)=step(d,D,y,-kdeg1*y*dt+gy2+gy3,dx,dt,N);
    u=U(:,i);
    v=V(:,i);
    w=W(:,i);
    y=Y(:,i);
    R1(:,i)=step2(u,r1,Rtot1,kon1,koff1,kdeg1,Rec,dt);
    R2(:,i)=step2(v,r2,Rtot2,kon2,koff2,kdeg2,Rec,dt);
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

function [n,flag]=grow(r1,r2,kg1,kg2,flag)
%this part is tricky. must return n integer and n<N
if flag==1,
    n=-1;
    return;
end
if sum(r1)<kg1 && sum(r2)>kg2,
    n=1;
elseif sum(r2)<kg2/2,
    n=-1;
    flag=1;
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


function u=step(d,D,u,F,dx,dt,N)
%I=eye(N,N);
I=eye(N-1,N-1);
%R=zeros(N,1);
u(1:N-1) =(I-d*(dt/dx^2)*D)\(u(1:N-1)+F(1:N-1));
%u =(I-d*(dt/dx^2)*D)\(u+F);
return

function [gv3,gy3]=rxn2(v,y,kr2,dt)
gv3=-v.*y*kr2*dt;
gy3=gv3;


function r=step2(u,r,rtot,kon,koff,kdeg,Rec,dt)
r=(kon*u(Rec).*(rtot-r)-(koff+kdeg)*r)*dt+r;
return
% function ru=step3(u,r,rtot,kon,koff,Rec,dt)
% ru=u*0;
% ru(Rec)=-(kon*u(Rec).*(rtot-r)-(koff))*dt;
% return