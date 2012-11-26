close all
U=zeros(1500,3);
U(1,:)=[0.01,0.01,0.01];
dt=0.001;
for i=2:15000,
    U(i,1)=(1-U(i-1,3))*U(i-1,1);
    U(i,2)=(1-U(i-1,1))*U(i-1,2);
    U(i,3)=(1-U(i-1,2))*U(i-1,3);
    U(i,:)=U(i,:)*dt+U(i-1,:);
end

t=dt:dt:dt*15000;


figure
plot(t,U(:,1),t,U(:,2),t,U(:,3))

