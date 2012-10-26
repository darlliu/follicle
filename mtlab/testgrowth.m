clc
flagin=1;
om=76;
on=6;
sm=71;
sn=1;
for i=1:160,
    [om,on,flagin]=growth2d(om,on,0,0,sm,sn,0,0,flagin);
    disp(om),disp(on),disp(flagin)
end
