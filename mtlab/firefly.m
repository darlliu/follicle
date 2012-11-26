% Fireflies sync to lamp
close all

n=100;
theta=random('normal',0,1,11,11);
[X,Y]=meshgrid(-10:2:10,-10:2:10);


Lpos=[0, 0];
L=0;

d=((X-Lpos(1)).^2+(Y-Lpos(2)).^2).^0.5;
%distance

fig1=figure(1);

surf(X,Y,theta);
axis([-10 10 -10 0 -1 1]);
view(0,90)
windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
Movie=moviein(40,fig1,windowsize);
for i=1:40,
    
    [L,theta]=step(theta, L, 0.8, d, 1, 0.5, @(x) x./5 );
    theta(6,6)=L;
    
    surf(X,Y,theta);
	axis([-10 10 -10 10 -0.5 0.5]);
    view(0,90)
    Movie(:,i)=getframe(fig1,windowsize);
end

movie(fig1, Movie, 30,3,windowsize);