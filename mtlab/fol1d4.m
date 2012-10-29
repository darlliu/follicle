function fol1d4
% predator prey model v1

clc
close all
N=200;
% discritize N pts along x
x=linspace (0,1,N);
dx=1/N;
% distance vector
% solve in [0,1]



M=8000;
% time points

dt=0.005; 
%discritized time;


r1=[0; 0];
r2=[0; 0];

R1=zeros(2,M);
R2=R1;
%receptor bound matrix
%initial at 0


Src=[10 11 12 13 14 15];
%the linear source
%not used in this model

Rec=[20 21];
%receptors source

Ends=[50 51];
%the initial ends;

%note: you can change the length of these vectors but remember to change R
%accordingly

u=zeros(N,1);
%prepare temp conc vectors
u(Ends)=1E-5;
%initial conditions, change to whatever.
v=u;

U=zeros(N,M);
U(:,1)=u;
V=U;

%total concentration-time matrix of ligands

D=MakeLaplacian1D(N);
% diffusion matrix

t=0:dt:dt*(M-1);
%time vector

%parameters:

Rtot1=1E-3;
Rtot2=Rtot1;
%max receptor amount
kon1=1;
koff1=1;
kon2=1;
koff2=1;
kdeg1=1E-3;
kdeg2=1E-3;
%rates
a1=1;
a2=1;
b1=1;
b2=1;
%coupling constants
d=1E-3;
%diffusion rate

for i= 2: M,
    gu=rxn(u,r1,kon1,koff1,Rtot1,Rec,dt,N);
    gv=rxn(v,r2,kon2,koff2,Rtot2,Rec,dt,N);
    Ends=Ends+grow(r1,r2);
    [gu2,gv2]=gen(r1,r2,a1,a2,b1,b2,Ends,N,dt);
    U(:,i)=step(d,D,u,gu+gu2,dx,dt,N);
    V(:,i)=step(d,D,v,gv+gv2,dx,dt,N);
    R1(:,i)=step2(u,r1,Rtot1,kon1,koff1,kdeg1,Rec);
    R2(:,i)=step2(v,r2,Rtot2,kon2,koff2,kdeg2,Rec);
    u=U(:,i);
    v=V(:,i);
    r1=R1(:,i);
    r2=R2(:,i);
end



%plotting routines
figure

plot(t,sum(U(Rec,:)),t,sum(V(Rec,:)));
title('receptor bound overtime')
legend('BMP_LR','WNT_LR')

fig1=figure;
plot(x,U(:,i),x,V(:,i))
legend('BMP','wnt')
windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
Movie=moviein(100,fig1,windowsize);
Movie(:,1)=getframe(fig1,windowsize);
frame=2;

for i=101:M/100:M,

    plot(x,U(:,i),x,V(:,i))
    legend('BMP','wnt')
    %plot(x,U(:,i))
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;

end

% size(Results)

% size(U)


% size(Stem)
 movie(fig1, Movie, 100,8,windowsize);
return

function n=grow(r1,r2)
%this part is tricky. must return n integer and n<N
n=0;
return
function gr=rxn(u,r,kon,koff,rtot,Rec,dt,N)
%gives reaction part given previous quantities and rates
gr=zeros(N,1);
gr(Rec)=dt*(-kon*u(Rec).*(rtot-r)+koff*(r));

return

function [gu2,gv2]=gen(r1,r2,a1,a2,b1,b2,Ends,N,dt)
gu2=zeros(N,1);
gv2=gu2;
gu2(Ends)=dt*(a1*(sum(r1))-a2*r1'*r2);
gv2(Ends)=dt*(b1*(r1'*r2)-b2*sum(r2));
return


function u=step(d,D,u,F,dx,dt,N)
I=eye(N-1,N-1);

%R=zeros(N,1);
u(1:N-1) =(I-d*(dt/dx^2)*D)\(u(1:N-1)+F(1:N-1));
return

function r=step2(u,r,rtot,kon,koff,kdeg,Rec)
r=kon*u(Rec).*(rtot-r)-(koff+kdeg)*r;
return
