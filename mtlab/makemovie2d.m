function Movie=makemovie2d (M,N,U,V,Is)
%movie
close all
[X, Y]=meshgrid(linspace(0,1,N), linspace(0,1,M));
%get grid
fig1=figure;
u=real(reshape(U(:,1), N,M)');
u=u./max(max(u));
v=real(reshape(V(:,1), N,M)');
v=v./max(max(v));
surf(X,Y,u,'FaceColor','interp','FaceLighting','phong');
hold on
surf(X,Y,v,'FaceColor','interp','FaceLighting','phong');
hold off
disp(max(max(u))/max(max(v)))
view(0,90)
windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
frame=2;
Movie=moviein(100,fig1,windowsize);
Movie(:,1)=getframe(fig1,windowsize);
for i=101:Is/100:Is,
    u=real(reshape(U(:,i), N,M)');
    u=u./max(max(u));
    v=real(reshape(V(:,i), N,M)');
    v=v./max(max(v));
    surf(X,Y,u,'FaceColor','interp','FaceLighting','phong');
    hold on
    surf(X,Y,v,'FaceColor','interp','FaceLighting','phong');
    hold off
    view(0,90)
    Movie(:,frame)=getframe(fig1,windowsize);
    frame=frame+1;
end
movie(fig1, Movie, 100,8,windowsize);
return