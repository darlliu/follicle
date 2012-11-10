function compare

% here we compare at least five numerical methods:
% 1. forward euler, 
% 2. explicit-implicit euler with a smaller explicit dt for nonlinear
% estimation
% 3. explicit-implicit with zero order explicit approximation (where
% f(n+1)~f(n))
% 4. rk4
% 5. edrk4 with dft space domain and rk4 time domain (above all central
% limit)


% our test function is x'=Dx+xy y'=Dy+xy with initial condition 1,1. 

clc
close all

x=linspace(0,1,20);
tmax=1;
% dt1=tmax/20;
% t1=linspace(0,tmax,20);
dt2=tmax/200;
t2=linspace(0,tmax,200);
dt3=tmax/2000;
t3=linspace(0,tmax,2000);
%our two dts, fwd will be very accurate!
N=20;
dx=1/20;
D=-2*eye(N,N)+diag(diag(eye(N-1,N-1)),1)+diag(diag(eye(N-1,N-1)),-1);
D(1,N)=1;
D(N,1)=1;
%periodic L

u=zeros(N,1);
v=u;
u(10)=1;
v(5)=1;
%our ic

%first fem
U=zeros(20,2000);
V=U;
U(:,1)=u;
V(:,1)=v;
for i=2:2000,
    U(:,i)=fem(D,U(:,i-1),V(:,i-1),dt3,dx);
    V(:,i)=fem(D,V(:,i-1),U(:,i-1),dt3,dx);

end
size(U)
U2=zeros(N,200);
V2=zeros(N,200);
E2=U2;
EE2=V2;
U2(:,1)=u;
V2(:,1)=v;

for i=2:200,
    U2(:,i)=bwd(D,U2(:,i-1),U(:,i*10),V(:,i*10),dt2,dx);
    V2(:,i)=bwd(D,V2(:,i-1),U(:,i*10),V(:,i*10),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end



figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);
size(X)
size(T)
h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Forward Euler, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Forward Euler,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);
size(X)
size(T)
h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Backward Euler, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Backward Euler,V');

subplot(3,2,5)
plot(t2,sum(E2))
title('Error, U');
subplot(3,2,6)
plot(t2,sum(EE2))

title('Error,V');

figure
plot(t3,U(3,:),t2,U2(3,:))
title('U3 plots')
legend('fwd','bwd')

%----------Here compares FEM with BWD--------%



for i=2:200,
    U2(:,i)=bwd2(D,U2(:,i-1),V2(:,i-1),dt2,dx);
    V2(:,i)=bwd2(D,V2(:,i-1),U2(:,i-1),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);
size(X)
size(T)
h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Forward Euler, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Forward Euler,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);
size(X)
size(T)
h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Backward Euler0, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Backward Euler0,V');

subplot(3,2,5)
plot(t2,sum(E2))
title('L2 sum Error, U');
subplot(3,2,6)
plot(t2,sum(EE2))
title('L2 sum Error,V');
xlabel('time')

figure
plot(t3,U(3,:),t2,U2(3,:))
title('U3 plots')
legend('fwd','bwd0')

%----------Here compares FEM with BWD0--------%
U3=zeros(size(U));
V3=U3;
U3(:,1)=u;
V3(:,1)=v;

for i=2:2000,
    [U3(:,i),V3(:,i)]=rk4(D,U3(:,i-1),V3(:,i-1),dt3,dx);

    E2(:,i)=abs(U(:,i)-U3(:,i));
    EE2(:,i)=abs(V(:,i)-V3(:,i));
end

figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);
size(X)
size(T)
h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Forward Euler, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Forward Euler,V');

subplot(3,2,3)
h=pcolor(T,X,U3);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('RK4, U');
subplot(3,2,4)
h=pcolor(T,X,V3);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('RK4,V');

subplot(3,2,5)
plot(t3,sum(E2))
title('L2 sum Error, U');
subplot(3,2,6)
plot(t3,sum(EE2))
title('L2 sum Error,V');
xlabel('time')

figure
plot(t3,U(3,:),t3,U3(3,:))
title('U3 plots')
legend('fwd','Rk4')

%----------Here compares FEM with RK4--------%
return


function unext=fem (D,u,v,dt,dx)
unext=((1/dx^2)*D*u+10*u.*v)*dt+u;
return

function unext=fem2(D,u,v,dt,dx)
unext=((1/dx^2)*D*u+10*u.*v)*dt;
return

function unext=bwd (D,u,u0,v0,dt,dx)
I=eye(20,20);
unext=(I-(dt/dx^2)*D)\(u+10*u0.*v0*dt);

function unext=bwd2 (D, u, v, dt,dx)
I=eye(20,20);
unext=(I-(dt/dx^2)*D)\(u+10*u.*v*dt);

function [unext,vnext]=rk4 (D,u,v,dt,dx)
u1=fem2(D,u,v,dt,dx);
v1=fem2(D,v,u,dt,dx);
%k1

u2=fem2(D,u+u1/2,v+v1/2,dt,dx);
v2=fem2(D,v+v1/2,u+u1/2,dt,dx);

%k2

u3=fem2(D,u+u2/2,v+v2/2,dt,dx);
v3=fem2(D,v+v2/2,u+u2/2,dt,dx);
%k3

u4=fem2(D,u+u3,v+v3,dt,dx);
v4=fem2(D,v+v3,u+u3,dt,dx);

%k4

unext=u+1/6*(u1+2*u2+2*u3+u4);
vnext=v+1/6*(v1+2*v2+2*v3+v4);
return
