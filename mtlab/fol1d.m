function fol1d
clc
close all
N=20;
x=linspace (0,1,N);
% solve in [0,1]
u=zeros(N,1);
u(1)=0.001;
v=u;
% discritize N pts along x
M=50000;
M2=M;
dt=0.001; %discritize time;
U=zeros(N,M);
V=U;
%result matrix
t=0:dt:dt*(M2-1);
%time vector
tnow=0;
for i= 2: M,
    [gu,gv]=gen(u,v,N);
    u=step(u,dt,N,gu);
    v=step(v,dt,N,gv);
    U(:,i)=u;
    V(:,i)=v;
    tnow=tnow+dt;
end


fig1=figure(1);
plot(x,U(:,1),x,V(:,1))
legend('BMP','wnt')
 windowsize=get(fig1,'Position');
 windowsize(1:2)=[0,0];
 Movie=moviein(100,fig1,windowsize);
 Movie(:,1)=getframe(fig1,windowsize);
frame=2;
for i=101:M/100:M,

    plot(x,U(:,i),x,V(:,i))
    legend('BMP','wnt')
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;
end
% size(Results)
% size(t)
% size(U)
% size(V)
% size(x)
% size(Stem)
 movie(fig1, Movie, 100,3,windowsize);



    

function [gu,gv]=gen (u,v,N)

gu=sum(u)*(1-1*sum(v)/N)/N;
gv=-sum(v)*(1-1*sum(u)/N)/N;

return


function unext=step (u, dt, N,generated)
%forward euler
dx=1/N;
D=-2*eye(N,N)+diag(diag(eye(N-1,N-1)),1)+diag(diag(eye(N-1,N-1)),-1);
D(1,1)=-1;
D(N-1,N-1)=-1;
D(N,N)=0;
D(N,N-1)=-1;
D(N,N-2)=1;
%D with mixed neumann condition

F=zeros(N,1);
F(1)=generated; 
F(N)=0;
%F(0) is the source of diffusion.


unext=(D*u+F*dx).*(dt/dx^2);

unext=unext+u;
return
