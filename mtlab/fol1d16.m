function fol1d16
%fol1d6 but with ifab2 method.
close all

global M N steps dt dx;
N=200;
M=4000;
dt=0.01;
steps=10;
% discritize N pts along x
x=linspace (0,1,N);
dx=1/N;
% distance vector
% solve in [0,1]



% time points

%discritized time;


r1=[0; 0];
r2=[0; 0];

R1=zeros(2,M);
R2=R1;
%receptor bound matrix
%initial at 0
global Src Rec Prc Ends; 

Src=[30 31];
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


t=0:dt:dt*(M-1);
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
kg1=0.8E-4;
kg2=1.5E-4;
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
D=MakeLaplacian1D(N);
Dm=expm((d*dt/dx^2)*D);
Dmm=expm((2*d*dt/dx^2)*D);
Dm2=expm((d2*dt/dx^2)*D);
Dmm2=expm((2*d2*dt/dx^2)*D);
Dm3=expm((d3*dt/dx^2)*D);
Dmm3=expm((2*d3*dt/dx^2)*D);
Dm4=expm((d4*dt/dx^2)*D);
Dmm4=expm((2*d4*dt/dx^2)*D);
% diffusion matrix
flag=0;
%rates
% a1=2;
% a2=2;
% b1=1;
% b2=1;
%coupling constants

H=zeros(M,1);
H(1)=50;
H(2)=50;
counter2=40;

[U(:,2),W(:,2),V(:,2),Y(:,2)]= fem(u,w,v,y,r1,r2,dt);
for i= 3: M,

    if counter==1,
        [n,flag]=grow(r1,r2,kg1,kg2,flag);
        Ends=Ends+n;
%         if Ends(2)>71,
%             Prc=Ends-2;
%         else
%             Prc=[];
%         end
        if Ends(2)>=90 && flag==1,
            Ends=[89,90];
            flag=1.5;
            counter2=40;
        end
        if flag==1.5,
            if counter2>1,
                Prc=Ends-2;
                counter2=counter2-1;
            else
                flag=2;
                Prc=[];
            end
        end
        
        if Ends(1)<=50,
            Ends=[50,51];
            flag=0;
        end
        counter=steps;
    end
    counter=counter-1;
    H(i)=Ends(1);
    
    
    [U(:,i),W(:,i),V(:,i),Y(:,i),R1(:,i),R2(:,i)]=...
        ifab2(U(:,i-1),W(:,i-1),V(:,i-1),Y(:,i-1),R1(:,i-1),R2(:,i-1), ...
        U(:,i-2),W(:,i-2),V(:,i-2),Y(:,i-2),R1(:,i-2),R2(:,i-2),dt);
%     u=U(:,i);
%     v=V(:,i);
%     w=W(:,i);
%     y=Y(:,i);
    r1=R1(:,i);
    r2=R2(:,i);
%     isnan(U(:,i))
%     isnan(v)
%     isnan(w)
%     isnan(y)
%     isnan(r1)
%     isnan(r2)
end



%plotting routines
figure
% subplot(1,2,1)
H2=H*max(max(R1))/50;
plot(t,sum(R1),t,sum(R2),t,H2);

title('receptor bound overtime')
legend('BMP_L_R','Wnt_L_R','Growth')
figure
subplot(3,1,1)
plot(H,sum(R1)/max(sum(R1)))
legend('BMP_s_s')
xlabel('Growth')
ylabel('Relative binding')
title('phase portrait')
subplot(3,1,2)
plot(H,sum(R2)/max(sum(R2)))
legend('Wnt_s_s')
xlabel('Growth')
ylabel('Relative binding')
title('phase portrait')
subplot(3,1,3)
plot(H,sum(R2)./(sum(R1)))
legend('Wnt/BMP_s_s')
xlabel('Growth')
ylabel('Relative binding')
% subplot(1,2,2)
% plot(t,sum(R2));
% title('receptor bound overtime')
% legend('Wnt_L_R')
% figure
% plot(t,H);
% title('growth vs time')
fig1=figure;
i=1;
% H(i);
% m1=max(max(U));
% m2=max(max(V));
% m3=max(max(W));
% m4=max(max(Y));
% m=max([m1,m2,m3,m4]);
%m=1E-3;
plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i),x(Rec),R1(:,i),'*',...
    x(Rec),R2(:,i),'o',...
    linspace(x(30),x(H(i)),10),linspace(0, U(H(i),i),10),'^k')
 %axis([0 1 0 m])
text(0,0,num2str(i*dt))
legend('BMP','wnt','Noggin','Dkk','BMP_L_R','WNT_L_R','Growth')
xlabel('Distance')
ylabel('Concentration')
windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
Movie=moviein(100,fig1,windowsize);
Movie(:,1)=getframe(fig1,windowsize);
frame=2;

for i=ceil(M/100)+1:M/100:M,

    plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i),x(Rec),R1(:,i),'*',...
        x(Rec),R2(:,i),'o',...
        linspace(x(30),x(H(i)),10),linspace(0, U(H(i),i),10),'^k')
    legend('BMP','wnt','Noggin','Dkk','BMP_L_R','WNT_L_R','Growth')
    xlabel('Distance')
    ylabel('Concentration')
    %axis([0 1 0 m])
    text(0,0,num2str(i*dt))
    %plot(x,U(:,i))
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;

end

% size(Results)

% size(U)


% size(Stem)
movie2avi(Movie,'fol1d15','fps',10)
 
return

function [n,flag]=grow(r1,r2,kg1,kg2,flag)
% sum(r1)
% sum(r2)
%this part is tricky. must return n integer and n<N
if flag==1,
    n=1;
    return;
elseif flag==2,
    n=-1;
elseif flag==1.5,
    n=0;
elseif sum(r1)<kg1 && sum(r2)>kg2,
    n=1;
    flag=1;
else
    n=0;
end

return

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
global N dx;

gu=rxn(u,r1,kon1,koff1,Rtot1,Rec,N);
gv=rxn(v,r2,kon2,koff2,Rtot2,Rec,N);
gu2=gen(kgen1*0.2,Src,N);
gv2=gen(kgen2*0.1,Ends,N);
gv2=gv2+gen(kgen2*0.1,Src,N);
gw2=gen(kgen1*0.6,Ends,N);
gy2=gen(kgen2*0.2,Src,N);

if isempty(Prc),
else
    gu2=gu2+gen(kgen1*0.5,Prc,N);
    gy2=gy2+gen(kgen2*0.9,Prc,N);
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
% function ru=step3(u,r,rtot,kon,koff,Rec,dt)
% ru=u*0;
% ru(Rec)=-(kon*u(Rec).*(rtot-r)-(koff))*dt;
% return