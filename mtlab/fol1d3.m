function fol1d3
clc
close all
N=20;
x=linspace (0,1,N);
% solve in [0,1]
u=zeros(N,1);
u(1)=0.0005;
v=u;
w=v;
y=w;
% discritize N pts along x
M=20000;
M2=M;
dt=0.002; %discritize time;
U=zeros(N,M);
V=U;
W=V;
Y=W;
%result matrix
t=0:dt:dt*(M2-1);
%time vector
n=7;
flag=0;
toplot=ones(1,M);
for i= 2: M,
    [gu,gv,gw,gy]=gen(u,v, w, y,n);
    u=step(u,w,dt,4,gu,N);
    w=step(w,u,dt,n,gw,N);
    v=step(v,y,dt,n,gv,N);
    y=step(y,v,dt,4,gy,N);
    U(:,i)=u;
    W(:,i)=w;
    V(:,i)=v;
    Y(:,i)=y;    
    toplot(i)=n;
    [n,flag]=nextn(u,v,n,flag);
end
size(U)
size(V)
size(t)
size(x)
figure
plot(toplot)
figure
plot(t,sum(U),t,sum(V));

legend('BMP','wnt')
fig1=figure;
plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
legend('BMP','wnt','Noggin','Dkk')
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
    plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
legend('BMP','wnt','Noggin','Dkk')
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
return

    

function [gu,gv,gw,gy]=gen (u,v,w,y,n)
gu=0.04*u(n)^2/(1+u(n)^2);

gw=0.05*w(n)^2/(1+w(n)^2);

gv=0.005*v(n)^2/(1+v(n)^2);

gy=0.01*y(n)^2/(1+y(n)^2);

return


function u=step(u,w,dt,n,gu,N)
%backward euler
dx=1/N;
D=-2*eye(N,N)+diag(diag(eye(N-1,N-1)),1)+diag(diag(eye(N-1,N-1)),-1);
D(1,1)=-1;
D(1,2)=0;
D(N-1,N-1)=-1;
D(N,N)=0;
D(N,N-1)=-1;
D(N,N-2)=1;
D(n,n)=-1;
D(n-1,n-1)=-1;
D(n,n+1)=1;
D(n-1,n-2)=1;
D=D*0.01;
%D with mixed neumann condition
F=zeros(N,1);
F(n)=gu/dx;
%F(0) is the source of diffusion.
I=eye(N,N);
R=-0.5*(u.*w./(1+u.*w)).*u-0.5*u;
u =(I-(dt/dx^2)*D)\(u+dt*F+R*dt);
%unext=(D*u+F).*(dt/dx^2);
return

        
