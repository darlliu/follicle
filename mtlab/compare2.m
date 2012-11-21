function compare2
clc
close all
global N M
N=300;
M=1000;
dt=0.01;
x=linspace (0,1,N);
t=0:dt:dt*(M-1);
err=[0 0 0 0 0 0];
el=[0 0 0 0 0 0];

[X,T]=meshgrid(t,x);
t=cputime;

U=fol1d13(M,N,dt);
%AB4BD4 is our reference: fourth order in both time and space
t2=cputime;
elapsed=t2-t;
t=t2;
figure
h=pcolor(T,X,U);
colormap(jet)
caxis([0 1E-6])
shading interp 
set(h,'edgecolor','none');
xlabel('position')
ylabel('time')
title('realplot, AB4BD4')
disp('AB4BD4 time')
disp('elapsed time:')
elapsed
el(1)=elapsed;

disp('Here are the sqrt(L2 error sum)')
E=0*U;

t=cputime;
V=fol1d9(M,N,dt);
t2=cputime;
elapsed=t2-t;
t=t2;
E=sqrt((U-V).^2);
figure
subplot(3,2,1)

h=pcolor(T,X,E);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 1E-6])
xlabel('position')
ylabel('time')
title('Error relative, AB2BD2');
disp('AB2BD2')
err(2)=sum(sum(E));
sum(sum(E))
disp('elapsed time:')
elapsed
el(2)=elapsed;
t=cputime;
V=fol1d8(M,N,dt);
t2=cputime;
elapsed=t2-t;
t=t2;


E=sqrt((U-V).^2);

subplot(3,2,2)

h=pcolor(T,X,E);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 1E-6])
xlabel('position')
ylabel('time')
title('Error relative, crn1');
disp('crn1')
err(3)=sum(sum(E));
sum(sum(E))
disp('elapsed time:')
elapsed
el(3)=elapsed;
t=cputime;
V=fol1d11(M,N,dt);
t2=cputime;
elapsed=t2-t;
t=t2;
E=sqrt((U-V).^2);

subplot(3,2,3)

h=pcolor(T,X,E);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 1E-6])
xlabel('position')
ylabel('time')
title('Error relative, ifab2');
disp('ifab2')
err(4)=sum(sum(E));
sum(sum(E))
disp('elapsed time:')
elapsed
el(4)=elapsed;
t=cputime;
V=fol1d12(M,N,dt);
t2=cputime;
elapsed=t2-t;
t=t2;
E=sqrt((U-V).^2);

subplot(3,2,4)

h=pcolor(T,X,E);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 1E-6])
xlabel('position')
ylabel('time')
title('Error relative, edrk2');
disp('edrk2')
err(5)=sum(sum(E));
sum(sum(E))
disp('elapsed time:')
elapsed
el(5)=elapsed;
t=cputime;
V=fol1d10(M,N,dt);
t2=cputime;
elapsed=t2-t;
t=t2;
E=sqrt((U-V).^2);

subplot(3,2,5)

h=pcolor(T,X,E);
colormap(jet)
shading interp 
set(h,'edgecolor','none');
caxis([0 1E-6])
xlabel('position')
ylabel('time')
title('Error relative, ifrk4');

disp('ifrk4')
el(6)=sum(sum(E));
sum(sum(E))
disp('elapsed time:')
elapsed
el(6)=elapsed;

figure
subplot(2,1,1)

bar(err)
title('error relative to AB4BD4')
ylim([-0.015 0.06])
text(1,-0.01,'AB4BD4')
text(2,-0.01,'AB2BD2')
text(3,-0.01,'CRN1')
text(4,-0.01,'IFAB2')
text(5,-0.01,'EDRK2')
text(6,-0.01,'IFRK4')
subplot (2,1,2)
bar(el)
title('solve time')
ylim([-10 120])
text(1,-5,'AB4BD4')
text(2,-5,'AB2BD2')
text(3,-5,'CRN1')
text(4,-5,'IFAB2')
text(5,-5,'EDRK2')
text(6,-5,'IFRK4')
