function Movie=makemovie2d (M,N,U,V,Is,t1,t2)
%movie

[X, Y]=meshgrid(linspace(0,1,N), linspace(0,1,M));
%get grid

fig1=figure;
u=real(reshape(U(:,1), N,M)');
v=real(reshape(V(:,1), N,M)');
% m1=max(max(u));
% u=u/m1;
% 
% m2=max(max(v));
% v=v/m2+1;
%right shift
subplot(1,2,1)
surf(X,Y,u,'FaceColor','interp','FaceLighting','phong','EdgeColor','None');
caxis([0 5E-4])
view(90,90)
title(t1)
xlabel('X')
ylabel('Y')
% shading interp
% colormap jet
colorbar
subplot(1,2,2)
surf(X,Y,v,'FaceColor','interp','FaceLighting','phong','EdgeColor','None');
caxis([0 5E-4])
view(90,90)
title(t2)
xlabel('X')
ylabel('Y')
% shading interp
% colormap jet
colorbar


windowsize=get(fig1,'Position');
windowsize(1:2)=[0,0];
frame=2;
Movie=moviein(100,fig1,windowsize);
Movie(:,1)=getframe(fig1,windowsize);
for i=101:Is/100:Is,
    u=real(reshape(U(:,i), N,M)');
    v=real(reshape(V(:,i), N,M)');
    
    subplot(1,2,1)
    surf(X,Y,u,'FaceColor','interp','FaceLighting','phong','EdgeColor','None');
    caxis([0 5E-4])
    view(90,90)
    xlabel('X')
    ylabel('Y')
    title(t1)
%     shading interp
%     colormap jet
    colorbar
    subplot(1,2,2)
    surf(X,Y,v,'FaceColor','interp','FaceLighting','phong','EdgeColor','None');
    caxis([0 5E-4])
    view(90,90)
    colorbar
    xlabel('X')
    ylabel('Y')
    title(t2)
%     shading interp
%     colormap jet



    Movie(:,frame)=getframe(fig1,windowsize);
    
    frame=frame+1;
end

return