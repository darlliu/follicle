%hw2 q1
close all

y=0:0.001:25;
r=(y.^2/(1+y.^2))./(y.*(1-y/20));
figure

plot(r,y,r,0*y, 'LineWidth',2)
axis([0 15 0 25])
title('bifurcation with positive r and Y');
r=linspace(0.4,2,5);
yp=zeros(size(y,2),5);
for i = 1:5,
    ypt=(y.*r(i).*(1-y/20))-(y.^2/(1+y.^2));
    yp(:,i)=ypt;
end
figure
plot(y,yp(:,1),y,yp(:,2),y,yp(:,3),y,yp(:,4),y,yp(:,5),y,zeros(size(y)),'LineWidth',2)
legend('r=0.4','r=0.8','r=1.2,','r=1.6','2','zero axis')
xlabel('population')
ylabel('population change')
title('steady states');
text(0.1,-13,'different r, three steady states, the one at 0 is not plotted for some reason.')
text(0.1,-14,'the second one from zero is unstable while the larger one is table')
axis([0 25 -15 15])
z=ones(100,2);
dt=0.05;
for i=2:100,
    y=z(i-1,:);
    z(i,:)=[y(1)*1.*(1-y(1)/20)-y(1).^2/(1+y(1).^2),y(2)*1.25*(1-y(2)/20)-y(2).^2/(1+y(2).^2)];
    z(i,:)=z(i-1,:)+z(i,:)*dt;
end
figure
t=dt:dt:dt*100;
size(t)
size(z)
plot(t',z(:,1),t',z(:,2), 'LineWidth',2)
title('population, r=1.25 compared to r=1');
legend('r=1','r=1.25')
xlabel('time')
ylabel('population')

figure
plot(t',z(:,2)-z(:,1), 'LineWidth',2)
title('increase of population from r=1.25 compared to r=1');
legend('increase')
xlabel('time')
ylabel('population')