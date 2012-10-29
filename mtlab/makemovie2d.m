function Movie=makemovie2d (M,N,U,Is)
%movie
close all
[X, Y]=meshgrid(linspace(0,1,N), linspace(0,1,M));
%get grid
fig1=figure;
u=reshape(U(:,1), N,M)';
surf(X,Y,u);
view(0,90)
windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
frame=2;
Movie=moviein(100,fig1,windowsize);
Movie(:,1)=getframe(fig1,windowsize);
for i=101:Is/100:Is,
    u=reshape(U(:,i), N,M)';
    surf(X,Y,abs(u));
    view(0,90)
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;
end
movie(fig1, Movie, 100,3,windowsize);
return