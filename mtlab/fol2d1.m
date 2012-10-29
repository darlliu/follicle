function [U,V,Y,W, Is]= fol2d1
clc
close all
M=20; %discritize y
N=15; %discritize x
dx=1/M/N;
% Note: even though M, N is set here, they cannot be changed
dt=0.005; % must use backward euler, dt should be ?
Is=10000; %discrete time pts
Src= [1+6*N, 1+11*N, 1+16*N]; %source 1 (stem cell level)
src=[1 6; 1 11; 1 16];
U=zeros(M*N, Is); %our dense reactant matrix
V=U;
W=V;
Y=W;
u=zeros(M*N,1); % the temp u
v=u;
w=v;
y=w;
I=eye(N*M, N*M);
O=zeros(N*M,1);
U(Src,1)=[0 0.001 0.0001]; % our initial condition
u(Src)=[0 0.001 0.0001];
t=0:dt:(Is-1)*dt;
%time vector
ends=[4 9; 4 14;  4 19];
flag=[0, 0, 0];
%growth flags

counter=500;
%growth countdown counters;
nabla=MakeLaplacian2D(N,M);
for i = 2:Is,
    counter=counter-1;
    if counter==1,
        [ends(1,2) ends(1,1) flag(1)]= ...
            growth2d(M,N,ends(1,2),ends(1,1),u, w, src(1,2), src(1,1), 2E-4, 1E-6, flag(1));
        [ends(2,2) ends(2,1) flag(2)]= ...
            growth2d(M,N,ends(2,2),ends(2,1),u, w, src(2,2), src(2,1), 2E-4, 1E-6, flag(2));
        [ends(3,2) ends(3,1) flag(3)]= ...
            growth2d(M,N,ends(3,2),ends(3,1),u, w, src(3,2), src(3,1), 2E-4, 1E-6, flag(3));  
        counter=500;
        ends
    end

    Ends=(ends(:,1)+(ends(:,2)-1)*N)';
    %growth part
    [gu,gv,gw,gy]=gen(u,v, w, y, Ends);
    %generation part
    U(:,i)=step (nabla, u , w, dt, dx, Src, gu, I, O, 0.002,0.4,0.1);
    W(:,i)=step (nabla, w , u, dt, dx, Ends, gw, I, O, 0.003,0.4,0.1);
    V(:,i)=step (nabla, v , y, dt, dx, Ends, gv, I, O, 0.0005,0.6,0.1);
    Y(:,i)=step (nabla, y , v, dt, dx, Src, gy, I, O, 0.002,0.6,0.1);
    %diffusion rxn part
    u=U(:,i);
    w=W(:,i);
    v=V(:,i);
    y=Y(:,i);    
end

% Plotting part
figure
% first, growth functions
% second, concentration vs t at shg
subplot(2,1,1)
Src=src(:,1)+(src(:,2)+3-1)*N;
plot(t,U(Src(1),:),t,U(Src(2),:),t,U(Src(3),:) ...
    );
legend('BMP, 1','BMP, 2','BMP,3 ');
subplot(2,1,2)
plot(t,V(Src(1),:),t,V(Src(2),:),t,V(Src(3),:) ...
    );
legend('Wnt, 1','Wnt, 2','Wnt,3 ');

title ('concentration at stem cell vs time')

return

function [gu,gv,gw,gy]=gen (u,v,w,y,n)
gu=0.05*u(n).^2./(1+u(n).^2)+0.005;

gw=0.02*w(n).^2./(1+w(n).^2)+0.004;

gv=0.005*v(n).^2./(1+v(n).^2)+0.001;

gy=0.01*y(n).^2./(1+y(n).^2)+0.001;

return


function u=step (nabla, u , v, dt, dx, Ends, gu, I, O, ds, rs ,des)
F=O;
D=nabla*ds*dx;

F(Ends)=gu;
R=-rs*sqrt((v).*u)-des*u;
u =(I-(dt/dx^2)*D)\(u+dt*F+R*dt);

return