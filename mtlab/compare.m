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
tmax=1.2;
dt1=tmax/1000;
t1=linspace(0,tmax,1000);
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



U3=zeros(size(U));
V3=U3;
E2=U3;
EE2=U3;
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

U=U3;
V=V3;
%RK4 is our best to beat

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

h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);

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
legend('rk4','bwd')

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

h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);

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
legend('rk4','bwd0')

%----------Here compares FEM with BWD0--------%




for i=2:200,
    U2(:,i)=bwd3(D,U2(:,i-1),V2(:,i-1),dt2,dx);
    V2(:,i)=bwd3(D,V2(:,i-1),U2(:,i-1),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);

h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);

h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Backward Euler1, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Backward Euler1,V');

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
legend('rk4','bwd1')

%-------here compares with bwd1------%
% 
% for i=2:200,
%     U2(:,i)=crn0(D,U2(:,i-1),V2(:,i-1),dt2,dx);
%     V2(:,i)=crn0(D,V2(:,i-1),U2(:,i-1),dt2,dx);
% 
%     E2(:,i)=abs(U(:,i*10)-U2(:,i)); EE2(:,i)=abs(V(:,i*10)-V2(:,i));
% end figure subplot(3,2,1) [X,T]=meshgrid(t3,x);
% 
% h=pcolor(T,X,U); colormap(jet) shading interp set(h,'edgecolor','none');
% caxis([0 0.14]) xlabel('position') ylabel('time') title('Runge Kutta,
% U'); subplot(3,2,2) h=pcolor(T,X,V); colormap(jet) shading interp
% xlabel('position') ylabel('time') set(h,'edgecolor','none'); caxis([0
% 0.14]) title('Runge Kutta,V');
% 
% subplot(3,2,3) [X,T]=meshgrid(t2,x);
% 
% h=pcolor(T,X,U2); colormap(jet) shading interp set(h,'edgecolor','none');
% caxis([0 0.14]) xlabel('position') ylabel('time') title('Cran-Nicolson0,
% U'); subplot(3,2,4) h=pcolor(T,X,V2); colormap(jet) shading interp
% xlabel('position') ylabel('time') set(h,'edgecolor','none'); caxis([0
% 0.14]) title('Crank-Nicolson0,V');
% 
% subplot(3,2,5) plot(t2,sum(E2)) title('L2 sum Error, U'); subplot(3,2,6)
% plot(t2,sum(EE2)) title('L2 sum Error,V'); xlabel('time')
% 
% figure plot(t3,U(3,:),t2,U2(3,:)) title('U3 plots') legend('rk4','crn0')

%----------Here compares FEM with CRN0--------%


for i=2:200,
    U2(:,i)=crn1(D,U2(:,i-1),V2(:,i-1),dt2,dx);
    V2(:,i)=crn1(D,V2(:,i-1),U2(:,i-1),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);

h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');

subplot(3,2,3)
[X,T]=meshgrid(t2,x);

h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Cran-Nicolson1, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Crank-Nicolson1,V');

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
legend('rk4','crn1')

%----------Here compares FEM with CRN1--------%

U2(:,1)=u;
V2(:,1)=v;
U2(:,2)=U(:,20);
V2(:,2)=V(:,20);

for i=3:200,
    U2(:,i)=abbd(D,U2(:,i-1),U2(:,i-2),V2(:,i-1),V2(:,i-2),dt2,dx);
    V2(:,i)=abbd(D,V2(:,i-1),V2(:,i-2),U2(:,i-1),U2(:,i-2),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);

h=pcolor(T,X,U);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');


subplot(3,2,3)
[X,T]=meshgrid(t2,x);

h=pcolor(T,X,U2);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('ABBD2, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('ABBD2,V');

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
legend('rk4','abbd2')

%----------Here compares FEM with ABBD2--------%

U2(:,1)=u;
V2(:,1)=v;
U2(:,2)=U(:,20);
V2(:,2)=V(:,20);

for i=3:200,
    U2(:,i)=ifab2(D,U2(:,i-1),U2(:,i-2),V2(:,i-1),V2(:,i-2),dt2,dx);
    V2(:,i)=ifab2(D,V2(:,i-1),V2(:,i-2),U2(:,i-1),U2(:,i-2),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);

h=pcolor(T,X,U);
colormap(jet)
%shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
%shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');


subplot(3,2,3)
[X,T]=meshgrid(t2,x);

h=pcolor(T,X,U2);
colormap(jet)
%shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('IFAB2, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
%shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('IFAB2,V');

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
legend('rk4','ifab2')

%----------Here compares FEM with IFAB2--------%



%explicit higher order method
U2=zeros(20,200);
V2=U2;
U2(:,1)=u;
V2(:,1)=v;


for i=2:200,
    [U2(:,i),V2(:,i)]=ifrk2(D,U2(:,i-1),V2(:,i-1),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);

h=pcolor(T,X,U);
colormap(jet)
%shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
%shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');


subplot(3,2,3)
[X,T]=meshgrid(t2,x);

h=pcolor(T,X,U2);
colormap(jet)
%shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('ifrk2, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
%shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('ifrk2,V');

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
legend('rk4','ifrk2')

%------IFRK2 revised------%

%explicit higher order method
U2=zeros(20,200);
V2=U2;
U2(:,1)=u;
V2(:,1)=v;
U2(:,2)=U(:,2);
V2(:,2)=V(:,2);

for i=2:200,
    [U2(:,i),V2(:,i)]=etd1(D,U2(:,i-1),V2(:,i-1),dt2,dx);

    E2(:,i)=abs(U(:,i*10)-U2(:,i));
    EE2(:,i)=abs(V(:,i*10)-V2(:,i));
end
figure
subplot(3,2,1)
[X,T]=meshgrid(t3,x);

h=pcolor(T,X,U);
colormap(jet)
%shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('Runge Kutta, U');
subplot(3,2,2)
h=pcolor(T,X,V);
colormap(jet)
%shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('Runge Kutta,V');


subplot(3,2,3)
[X,T]=meshgrid(t2,x);

h=pcolor(T,X,U2);
colormap(jet)
%shading interp 
set(h,'edgecolor','none');
caxis([0 0.14])
xlabel('position')
ylabel('time')
title('etd1, U');
subplot(3,2,4)
h=pcolor(T,X,V2);
colormap(jet)
%shading interp 
xlabel('position')
ylabel('time')
set(h,'edgecolor','none');
caxis([0 0.14])
title('etd1,V');

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
legend('rk4','etd1')

%----------Here compares FEM with IFRK2--------%


% %close all U2=zeros(20,1000); V2=U2; U2(:,1)=u; V2(:,1)=v;
% U2(:,2)=U3(:,20); V2(:,2)=V3(:,20); E2=U2; EE2=V2; for i=3:1000, %
% U2(:,i)=sbdf(D,U2(:,i-1),U2(:,i-2),... %
% V2(:,i-1),V2(:,i-2),dt2,dx); % %
% V2(:,i)=sbdf(D,V2(:,i-1),V2(:,i-2),... %
% U2(:,i-1),U2(:,i-2),dt2,dx);
%     U2(:,i)=sbdf(D,U2(:,i-1),U2(:,i-2),...
%         V2(:,i-1),dt1,dx);
% 
%     V2(:,i)=sbdf(D,V2(:,i-1),V2(:,i-2),...
%         U2(:,i-1),dt1,dx);
% 
%     E2(:,i)=abs(U(:,i*2)-U2(:,i)); EE2(:,i)=abs(V(:,i*2)-V2(:,i));
% end figure subplot(3,2,1) [X,T]=meshgrid(t3,x);
% 
% h=pcolor(T,X,U); colormap(jet) shading interp set(h,'edgecolor','none');
% caxis([0 0.14]) xlabel('position') ylabel('time') title('Runge Kutta,
% U'); subplot(3,2,2) h=pcolor(T,X,V); colormap(jet) shading interp
% xlabel('position') ylabel('time') set(h,'edgecolor','none'); caxis([0
% 0.14]) title('Runge Kutta,V');
% 
% subplot(3,2,3) [X,T]=meshgrid(t1,x);
% 
% h=pcolor(T,X,U2); colormap(jet) shading interp set(h,'edgecolor','none');
% caxis([0 0.14]) xlabel('position') ylabel('time') title('SBDF2 , U');
% subplot(3,2,4) h=pcolor(T,X,V2); colormap(jet) shading interp
% xlabel('position') ylabel('time') set(h,'edgecolor','none'); caxis([0
% 0.14]) title('SBDF2 ,V');
% 
% 
% 
% subplot(3,2,5) plot(t1,sum(E2)) title('L2 sum Error, U'); subplot(3,2,6)
% plot(t1,sum(EE2)) title('L2 sum Error,V'); xlabel('time')
% 
% figure plot(t3,U(3,:),t1,U2(3,:)) title('U3 plots') legend('rk4','SBDF2')
% 
% %----------Here compares FEM with SBDF2--------%

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

function unext=bwd3 (D, u, v, dt,dx)
%implicit-explicit euler, with 1st order explicit approximation
I=eye(20,20);
u2=fem(D,u,v,dt,dx);
v2=fem(D,v,u,dt,dx);
unext=(I-(dt/dx^2)*D)\(u+10*u2.*v2*dt);

function unext=crn0 (D, u, v, dt,dx)
%implicit-explicit euler, with 1st order explicit approximation
I=eye(20,20);
unext=(I-(dt/dx^2)*D)\(u+0.5*dt*(2*10*u.*v+(1/dx^2)*(D*u)));

function unext=crn1 (D, u, v, dt,dx)
%implicit-explicit euler, with 1st order explicit approximation
I=eye(20,20);
u2=fem(D,u,v,dt,dx);
v2=fem(D,v,u,dt,dx);
unext=(I-(dt/dx^2)*D)\(u+0.5*dt*(10*u2.*v2+10*u.*v+(1/dx^2)*(D*u)));

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
% 
% function unext=abbd2 (D,u,u1,v,v1,dt,dx)
% %AB4BD4 IMEX.
% %the first 3 entries are loaded from RK4.
% unext=(25-12*(dt/dx^2)*D)\(48*u-36*u1+16*u2-3*u3+48*dt*u.*v...
%     -72*dt*u1.*v1+48*dt*u2.*v2-12*dt*u3.*v3);
% 
% function unext=sbdf(D,u,u1,v,v1,dt,dx)
% unext = (3-2*(dt/dx^2)*D)\(2*dt*(2*10*u.*v-10*u1.*v1)+4*u-u1);
% return

function unext=sbdf(D,u,u1,v,dt,dx)
I=eye(20,20);
unext = (I-(dt/dx^2)*D)\(2*dt*(10*u.*v)+(dt/dx^2)*D*u1+u1);
return

function unext=abbd(D,u,u1,v,v1,dt,dx)
I=eye(20,20);
unext = (3*I-2*(dt/dx^2)*D)\(4*u-u1+4*dt*(10*u.*v)-2*dt*(10*u1.*v1));
return

function unext=ifab2(D,u,u1,v,v1,dt,dx)

Dm=expm((dt/dx^2)*D);
Dm2=expm((2*dt/dx^2)*D);
unext=Dm*u+3*dt/2*Dm*(10*u.*v)-dt/2*Dm2*(10*u1.*v1);

function [unext,vnext]=ifrk2(D,u,v,dt,dx)
Dm=expm((dt/dx^2)*D);
u1=Dm*(u+dt*10*u.*v);
v1=Dm*(v+dt*10*u.*v);
f1=10*u1.*v1;
unext=Dm*u+dt/2*(Dm*(10*u.*v)+f1);
vnext=Dm*v+dt/2*(Dm*(10*u.*v)+f1);

function [unext,vnext]=etd2rk(D,u,v,u1,v1,dt,dx)
I=eye(20,20);
%svd
[lu,ls,lv]=svd((1/dx^2)*D);
ls(find(ls~=0))=1./ls(find(ls~=0));
%singular inverse of s
L=inv(lv)*ls*inv(lu);
Dm=expm((dt/dx^2)*D);
M1=L*(Dm-I)+dt*Dm*(I-L*((1/dx^2)*D));
M2=L*L*(Dm-(I+(dt/(dx^2))*D))+(1/2)*dt^2*Dm*(I-L*((1/dx^2)*D));
ua=Dm*u+M1*(10*u.*v);
va=Dm*v+M1*(10*u.*v);
unext=ua+M2*(10*u.*v-10*u1.*v1)/dt;
vnext=va+M2*(10*u.*v-10*u1.*v1)/dt;



function [unext,vnext]=etd1(D,u,v,dt,dx)
I=eye(20,20);
%svd
[lu,ls,lv]=svd((1/dx^2)*D);
ls(find(ls~=0))=1./ls(find(ls~=0));
%singular inverse of s
L=inv(lv)*ls*inv(lu);
Dm=expm((dt/dx^2)*D);
M1=L*(Dm-I)+dt*Dm*(I-L*((1/dx^2)*D));
%M2=L*L*(Dm-(I-dt*L))+(1/2)*dt^2*Dm*(I-L*((1/dx^2)*D));
%ua=Dm*u+M1*(10*u.*v);
%va=Dm*v+M1*(10*u.*v);
%unext=ua+M2*(10*ua.*va-10*u.*v)/dt;
%vnext=va+M2*(10*ua.*va-10*u.*v)/dt;

unext=Dm*u+M1*(10*u.*v);
vnext=Dm*v+M1*(10*u.*v);