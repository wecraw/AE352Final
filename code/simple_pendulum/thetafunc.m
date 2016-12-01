function  out  = thetafunc ( time, pos )

l=30;%pendulum length in meters
g=9.8;%m/s^2
theta=pos(1);
thetadot=pos(2);

thetadoubledot=-g*sin(theta)/l;
 
out=[thetadot; thetadoubledot];
end
