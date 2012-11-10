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
dt1=tmax/20;
t1=linspace(0,tmax,20);
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
D=D;

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
title('Forward Euler,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);
size(X)
size(T)
h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
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
title('Backward Euler,V');
sum(E2)
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
title('Forward Euler,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);
size(X)
size(T)
h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
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
title('Backward Euler0,V');
sum(E2)
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

return


function unext=fem (D,u,v,dt,dx)
unext=((1/dx^2)*D*u+10*u.*v)*dt+u;
return

function unext=bwd (D,u,u0,v0,dt,dx)
I=eye(20,20);
unext=(I-(dt/dx^2)*D)\(u+10*u0.*v0*dt);

function unext=bwd2 (D, u, v, dt,dx)
I=eye(20,20);
unext=(I-(dt/dx^2)*D)\(u+10*u.*v*dt);