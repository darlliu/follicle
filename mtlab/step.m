function [L,theta]=step (theta0, L0, dt, D, w, omega, fx)

L=L0+dt*omega;

theta=(w+fx(D).*sin(L-theta0))*dt;


