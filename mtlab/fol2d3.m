function [U,V,Y,W, R1, R2, Is2]=fol2d3
%fol1d6 but with ifab2 method.
clc
clear all
close all
global dx dy;

% discritize N pts along 
global Is N M steps dt;

% distance vector
% solve in [0,1]

Is=4500;
M=20; %discritize y
N=20; %discritize x
dx=1/N;
dy=1/M;
steps=20;
% time points unfortunately now have to set M==N
dt=0.01;
%discritized time;


r1=[0; 0; 0];
r2=r1;
%three receptors

R1=zeros(3,M);
R2=R1;
%receptor bound matrix
%initial at 0
global Src src Rec rec Prc prc Ends ends; 

Src= [1+6*N, 1+11*N, 1+16*N]; %source 1 (stem cell level)
src=[1 6; 1 11; 1 16];
%the linear source

Rec= [3+8*N, 3+13*N, 3+18*N];
rec=[3, 8; 3, 13; 3, 18];
%receptors source

Prc=[];
prc=[];
%precortex

ends=[4 9; 4 14;  4 19];
Ends=(ends(:,1)+(ends(:,2)-1)*N)';
%the initial ends;

%note: you can change the length of these vectors but remember to change R
%accordingly

u=zeros(M*N,1);
%prepare temp conc vectors
u(Ends)=1E-5;
%initial conditions, change to whatever.
v=u;
w=v;
y=w;
U=zeros(M*N, Is); %our dense reactant matrix
U(:,1)=u; 
V=U;
W=V;
Y=W;
%total concentration-time matrix of ligands


t=0:dt:dt*(Is-1);
%time vector

%parameters:
global d d2 d3 d4 D Dm Dmm Dm2 Dmm2 Dm3 Dmm3 Dm4 Dmm4 ...
    Rtot1 Rtot2 kon1 koff1 kon2 koff2 kdeg1 kdeg2 kg1 kg2...
    kgen1 kgen2 kr1 kr2 counter;

Rtot1=1E-1;
Rtot2=Rtot1;
%max receptor amount
kon1=.5;
koff1=.2;
kon2=.5;
koff2=.2;
kdeg1=1E-1;
kdeg2=1E-1;
kg1=11E-5;
kg2=2.3E-5;
kgen1=2.1E-4;
kgen2=2E-4;
kr1=5;
kr2=5;
d=5E-2;
d2=1E-3;
d3=2E-1;
d4=1E-1;
%diffusion rate

counter=steps;
D=MakeLaplacian2D(N,M);
Dm=expm((d*dt/dx^2)*D);
Dmm=expm((2*d*dt/dx^2)*D);
Dm2=expm((d2*dt/dx^2)*D);
Dmm2=expm((2*d2*dt/dx^2)*D);
Dm3=expm((d3*dt/dx^2)*D);
Dmm3=expm((2*d3*dt/dx^2)*D);
Dm4=expm((d4*dt/dx^2)*D);
Dmm4=expm((2*d4*dt/dx^2)*D);
% diffusion matrix
flag=[0, 0, 0];
%rates
% a1=2;
% a2=2;
% b1=1;
% b2=1;
%coupling constants

H=zeros(M,1);
[U(:,2),W(:,2),V(:,2),Y(:,2)]= fem(u,w,v,y,r1,r2,dt);
for i= 3: Is,
    counter=counter-1;
    if counter==1,
        [ends(1,2) ends(1,1) flag(1)]= ...
            growth2d2(M,N,ends(1,2),ends(1,1),r1(1), r2(1), src(1,2), src(1,1), kg1, kg2, flag(1));
        [ends(2,2) ends(2,1) flag(2)]= ...
            growth2d2(M,N,ends(2,2),ends(2,1),r1(2), r2(2), src(2,2), src(2,1), kg1, kg2, flag(2));
        [ends(3,2) ends(3,1) flag(3)]= ...
            growth2d2(M,N,ends(3,2),ends(3,1),r1(3), r2(3), src(3,2), src(3,1), kg1, kg2, flag(3));  
        counter=steps;
        ends
        Ends=(ends(:,1)+(ends(:,2)-1)*N)';
        if flag(1)~=0,
            prc(1,:)=ends(1,:)-1;
            Prc=abs(prc(:,1)+(prc(:,2)-1)*N)';
        end
        if flag(2)~=0,
            prc(2,:)=ends(2,:)-1;
            Prc=abs(prc(:,1)+(prc(:,2)-1)*N)';
        end
        if flag(3)~=0,
            prc(3,:)=ends(3,:)-1;
            Prc=abs(prc(:,1)+(prc(:,2)-1)*N)';
        end

        Src
        Ends
        Rec
        Prc
    end
    
    [U(:,i),W(:,i),V(:,i),Y(:,i),R1(:,i),R2(:,i)]=ifab2(u,w,v,y,r1,r2, ...
    U(:,i-2),W(:,i-2),V(:,i-2),Y(:,i-2),R1(:,i-2),R2(:,i-2),dt);
    u=U(:,i);
    v=V(:,i);
    w=W(:,i);
    y=Y(:,i);
    r1=R1(:,i);
    r2=R2(:,i);

end


% Plotting part
figure
% first, growth functions
% second, concentration vs t at shg
subplot(2,1,1)
Src=src(:,1)+(src(:,2)+3-1)*N;
plot(t,R1(1,:),t,R1(2,:),t,R1(3,:) ...
    );
legend('BMP, 1','BMP, 2','BMP,3 ');
subplot(2,1,2)
plot(t,R2(1,:),t,R2(2,:),t,R2(3,:) ...
    );
legend('Wnt, 1','Wnt, 2','Wnt,3 ');

title ('concentration at stem cell vs time')

% for i=101:M/100:M,
% 
%     plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
%     legend('BMP','wnt','Noggin','Dkk')
%     %plot(x,U(:,i))
%     Movie(:,frame)=getframe(fig1,windowsize);
%     frame=frame+1;
% 
% end

% size(Results)

% size(U)


% % size(Stem)
%  movie(fig1, Movie, 100,8,windowsize);
 Is2=Is;
return
% 
% function [n,flag]=grow(r1,r2,kg1,kg2,flag)
% %this part is tricky. must return n integer and n<N
% if flag==1,
%     n=1;
%     return;
% elseif flag==2,
%     n=-1;
% elseif sum(r1)<kg1 && sum(r2)>kg2,
%     n=1;
%     flag=1;
% else
%     n=0;
% end
% 
% return

function [un,wn,vn,yn,r1n,r2n]=ifab2(u,w,v,y,r1,r2,u2,w2,v2,y2,r12,r22,dt)
global d d2 d3 d4 D Dm Dmm Dm2 Dmm2 Dm3 Dmm3 Dm4 Dmm4
[uf,wf,vf,yf,r1f,r2f]=bigf(u,w,v,y,r1,r2,dt);
[uf2,wf2,vf2,yf2,r1f2,r2f2]=bigf(u2,w2,v2,y2,r12,r22,dt);

un=Dm*u+3/2*Dm*uf-1/2*Dmm*uf2;
wn=Dm3*w+3/2*Dm3*wf-1/2*Dmm3*wf2;
vn=Dm2*v+3/2*Dm2*vf-1/2*Dmm2*vf2;
yn=Dm4*y+3/2*Dm4*yf-1/2*Dmm4*yf2;
r1n=r1+3/2*r1f-1/2*r1f2;
r2n=r2+3/2*r2f-1/2*r2f2;
return



function [u2,w2,v2,y2,r12,r22]=fem(u,w,v,y,r1,r2,dt)
%used as an approximation
global D d d2 d3 d4 dx;
[uf,wf,vf,yf,r1f,r2f]=bigf(u,w,v,y,r1,r2,dt);

u2=(dt/dx^2)*d*D*u+uf*dt+u;
v2=(dt/dx^2)*d2*D*v+vf*dt+v;
w2=(dt/dx^2)*d3*D*w+wf*dt+w;
y2=(dt/dx^2)*d4*D*y+yf*dt+y;
r12=r1+r1f*dt;
r22=r2+r2f*dt;
return


function [uf,wf,vf,yf,r1f,r2f]=bigf(u,w,v,y,r1,r2,dt)
global d D Rtot1 Rtot2 kon1 koff1 kon2 koff2 kdeg1 kdeg2 kg1 kg2...
    kgen1 kgen2 kr1 kr2 counter;
global Src Rec Prc Ends; 
global N M dx;

gu=rxn(u,r1,kon1,koff1,Rtot1,Rec,N*M);
gv=rxn(v,r2,kon2,koff2,Rtot2,Rec,N*M);
gu2=gen(kgen1*0.2,Src,N*M);
gv2=gen(kgen2*0.1,Ends,N*M);
gv2=gv2+gen(kgen2*0.1,Src,N*M);
gw2=gen(kgen1*0.6,Ends,N*M);
gy2=gen(kgen2*0.2,Src,N*M);

if isempty(Prc),
else
    gu2=gu2+gen(kgen1*0.5,Prc,N*M);
    gy2=gy2+gen(kgen2*0.9,Prc,N*M);
end
[gu3,gw3]=rxn2(u,w,kr1);
[gv3,gy3]=rxn2(v,y,kr2);
%  111
% mean(u)
% mean(abs(gu))
% mean(gu2)
% mean(abs(gu3))
% mean(kdeg1*u)
uf=dt*(-kdeg1*u+gu+gu2+gu3);
vf=dt*(-kdeg2*v+gv+gv2+gv3);
wf=dt*(-kdeg1*w+gw2+gw3);
yf=dt*(-kdeg1*y+gy2+gy3);
r1f=step2(u,r1,Rtot1,kon1,koff1,kdeg1,Rec,dt);
r2f=step2(v,r2,Rtot2,kon2,koff2,kdeg2,Rec,dt);
return

function gr=rxn(u,r,kon,koff,rtot,Rec,N)
global dt
%gives reaction part given previous quantities and rates
gr=zeros(N,1);
gr(Rec)=((-kon/dt)*u(Rec).*(rtot-r)+(koff/dt)*(r));

return

function gu2=gen(kgen1,Src,N)
global dt
gu2=zeros(N,1);
gu2(Src)=kgen1/dt;

%gu2(Src)=dt*(sum(r1)/(kgen1+sum(r1)));
%gv2(Ends)=dt*(sum(r2)/(kgen2+sum(r2)));
return

% 
% function u=step(d,D,u,F,dx,dt,N)
% %I=eye(N,N);
% I=eye(N-1,N-1);
% %R=zeros(N,1);
% u(1:N-1) =(I-d*(dt/dx^2)*D)\(u(1:N-1)+F(1:N-1));
% %u =(I-d*(dt/dx^2)*D)\(u+F);
% return

function [gv3,gy3]=rxn2(v,y,kr2)
global dt
gv3=-v.*y*kr2/dt;
gy3=gv3;


function r=step2(u,r,rtot,kon,koff,kdeg,Rec,dt)

r=((kon/dt)*u(Rec).*(rtot-r)-(koff/dt+kdeg/dt)*r)*dt;
return