function fol1d3
clc
close all
N=20;
x=linspace (0,1,N);
% solve in [0,1]
u=zeros(N,1);
u(1)=0.001;
v=u;
% discritize N pts along x
M=10000;
M2=M;
dt=0.01; %discritize time;
U=zeros(N,M);
V=U;
%result matrix
t=0:dt:dt*(M2-1);
%time vector
tnow=0;
n=1;
toplot=ones(1,M);
for i= 2: M,
    [gu,gv]=gen(u,v, n);
    u=step(u,dt,N,gu);
    v=step(v,dt,N,gv);
    U(:,i)=u;
    V(:,i)=v;
    n=next(u,v,n);
    toplot(i)=n;
end


fig1=figure(1);
plot(x,U(:,1),x,V(:,1))
legend('BMP','wnt')
 windowsize=get(fig1,'Position');
 windowsize(1:2)=[0,0];
 Movie=moviein(100,fig1,windowsize);
 Movie(:,1)=getframe(fig1,windowsize);
frame=2;
results=zeros(2,100);
results(:,1)=[0.001 ; 0.001];
j=2;
for i=101:M/100:M,
    results(1,j)=u(toplot(j));
    results(2,j)=v(toplot(j));
    plot(x,U(:,i),x,V(:,i))
    legend('BMP','wnt')
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;
    j=j+1;
end

% size(Results)
 size(t)
% size(U)
 size(toplot)
 size(results)
% size(Stem)
 movie(fig1, Movie, 100,3,windowsize);

    

function [gu,gv]=gen (u,v,n)
gu=u(n)*(1-1*v(n));
gv=-v(n)*(1-1*u(n));
return

function nnext=next(u,v,nprev)
nnext=ceil((v(nprev)-u(nprev))+nprev)
if nnext<=0,
    nnext=1;
end
if nnext>20,
    nnext=20;
end
return

function unext=step (u, dt, N,generated)
%backward euler
dx=1/N;
D=-2*eye(N,N)+diag(diag(eye(N-1,N-1)),1)+diag(diag(eye(N-1,N-1)),-1);
D(1,1)=-1;
D(N-1,N-1)=-1;
D(N,N)=0;
D(N,N-1)=-1;
D(N,N-2)=1;
%D with mixed neumann condition

F=zeros(N,1);
F(1)=generated/dx; 
%F(0) is the source of diffusion.
I=eye(N,N);
unext =(I-(dt/dx^2)*D)\(u+ dt*F);
%unext=(D*u+F).*(dt/dx^2);

return
