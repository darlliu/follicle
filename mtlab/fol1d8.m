function U=fol1d8(M,N,dt)
%fol1d6 but with crn1



global  dx;

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
global d D Rtot1 Rtot2 kon1 koff1 kon2 koff2 kdeg1 kdeg2 kg1 kg2...
    kgen1 kgen2 kr1 kr2 counter;

Rtot1=1E-3;
Rtot2=Rtot1;
%max receptor amount
kon1=5;
koff1=2;
kon2=5;
koff2=2;
kdeg1=1E-3;
kdeg2=1E-3;
kg1=10;
kg2=5E-10;
kgen1=2.1E-5;
kgen2=2E-5;
kr1=14;
kr2=140;
d=1E-1;
%diffusion rate

counter=20;
D=MakeLaplacian1D(N);
% diffusion matrix
flag=0;
%rates
% a1=2;
% a2=2;
% b1=1;
% b2=1;
%coupling constants

H=zeros(M,1);
for i= 2: M,

    if counter==1,
        [n,flag]=grow(r1,r2,kg1,kg2,flag);
        Ends=Ends+n;
        if Ends(2)>65,
            Prc=Ends-2;
        else
            Prc=[];
        end
        if Ends(2)>=90,
            Ends=[89,90];
            flag=2;
        end
        if Ends(1)<=50,
            Ends=[50,51];
            flag=0;
        end
        counter=20;
    end
    counter=counter-1;
    H(i)=Ends(1);


    
    
    [U(:,i),W(:,i),V(:,i),Y(:,i),R1(:,i),R2(:,i)]=crn1(u,w,v,y,r1,r2,dt);
    u=U(:,i);
    v=V(:,i);
    w=W(:,i);
    y=Y(:,i);
    r1=R1(:,i);
    r2=R2(:,i);

end

% 
% 
% %plotting routines
% figure
% subplot(1,2,1)
% plot(t,sum(R1));
% title('receptor bound overtime')
% legend('BMP_L_R')
% subplot(1,2,2)
% plot(t,sum(R2));
% title('receptor bound overtime')
% legend('Wnt_L_R')
% figure
% plot(t,H);
% title('growth vs time')
% fig1=figure;
% plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
% legend('BMP','wnt','Noggin','Dkk')
% windowsize=get(fig1,'Position');
% windowsize(1:2)=[0,0];
% Movie=moviein(100,fig1,windowsize);
% Movie(:,1)=getframe(fig1,windowsize);
% frame=2;
% 
% for i=101:M/100:M,
% 
%     plot(x,U(:,i),x,V(:,i),x,W(:,i),x,Y(:,i))
%     legend('BMP','wnt','Noggin','Dkk')
%     %plot(x,U(:,i))
%     Movie(:,frame)=getframe(fig1,windowsize);
%     frame=frame+1;
% 
% end
% 
% % size(Results)
% 
% % size(U)
% 
% 
% % size(Stem)
%  movie(fig1, Movie, 100,8,windowsize);
 
return

function [n,flag]=grow(r1,r2,kg1,kg2,flag)
%this part is tricky. must return n integer and n<N
if flag==1,
    n=1;
    return;
elseif flag==2,
    n=-1;
elseif sum(r1)<kg1 && sum(r2)>kg2,
    n=1;
    flag=1;
else
    n=0;
end

return

function [un,wn,vn,yn,r1n,r2n]=crn1(u,w,v,y,r1,r2,dt)
global d D dx N;
I=eye(N,N);
[uf,wf,vf,yf,r1f,r2f]=bigf(u,w,v,y,r1,r2,dt);
[u2,w2,v2,y2,r12,r22]=fem(u,w,v,y,r1,r2,dt);
[uf2,wf2,vf2,yf2,r1f2,r2f2]=bigf(u2,w2,v2,y2,r12,r22,dt);
un=(I-d*(dt/dx^2)*D)\(u+0.5*(uf+uf2+(dt/dx^2)*(d*D*u)));
vn=(I-d*(dt/dx^2)*D)\(v+0.5*(vf+vf2+(dt/dx^2)*(d*D*v)));
wn=(I-d*(dt/dx^2)*D)\(w+0.5*(wf+wf2+(dt/dx^2)*(d*D*w)));
yn=(I-d*(dt/dx^2)*D)\(y+0.5*(yf+yf2+(dt/dx^2)*(d*D*y)));
r1n=r1+r1f;
r2n=r2+r2f;
return

function [u2,w2,v2,y2,r12,r22]=fem(u,w,v,y,r1,r2,dt)
%used as an approximation
global D d dx;
[uf,wf,vf,yf,r1f,r2f]=bigf(u,w,v,y,r1,r2,dt);

u2=(dt/dx^2)*d*D*u+uf*dt+u;
v2=(dt/dx^2)*d*D*v+vf*dt+v;
w2=(dt/dx^2)*d*D*w+wf*dt+w;
y2=(dt/dx^2)*d*D*y+yf*dt+y;
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
gv2=gen(kgen2*0.3,Ends,N);
gv2=gv2+gen(kgen2*0.3,Src,N);
gw2=gen(kgen1*0.9,Ends,N);
gy2=gen(kgen2*0.2,Src,N);
if isempty(Prc),
else
    gu2=gu2+gen(kgen1*0.5,Prc,N);
    gy2=gy2+gen(kgen2*0.9,Prc,N);
end
[gu3,gw3]=rxn2(u,w,kr1);
[gv3,gy3]=rxn2(v,y,kr2);
uf=dt*(-kdeg1*u+gu+gu2+gu3);
vf=dt*(-kdeg2*v+gv+gv2+gv3);
wf=dt*(-kdeg1*w+gw2+gw3);
yf=dt*(-kdeg1*y+gy2+gy3);
r1f=step2(u,r1,Rtot1,kon1,koff1,kdeg1,Rec,dt);
r2f=step2(v,r2,Rtot2,kon2,koff2,kdeg2,Rec,dt);
return

function gr=rxn(u,r,kon,koff,rtot,Rec,N)
%gives reaction part given previous quantities and rates
gr=zeros(N,1);
gr(Rec)=(-kon*u(Rec).*(rtot-r)+koff*(r));

return

function gu2=gen(kgen1,Src,N)
gu2=zeros(N,1);
gu2(Src)=kgen1;

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
gv3=-v.*y*kr2;
gy3=gv3;


function r=step2(u,r,rtot,kon,koff,kdeg,Rec,dt)
r=(kon*u(Rec).*(rtot-r)-(koff+kdeg)*r)*dt;
return
% function ru=step3(u,r,rtot,kon,koff,Rec,dt)
% ru=u*0;
% ru(Rec)=-(kon*u(Rec).*(rtot-r)-(koff))*dt;
% return