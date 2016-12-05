function  out  = thetafunc ( time, pos )
g=9.8;
l=20;
r=12;
theta=pos(1);
thetadot=pos(2);
thetadoubledot=-g*sin(theta)/l;
phi=pos(3);
phidot=pos(4);
phidoubledot=0;

%Rotation Matrix Stuff
% O_1to2=[0;l*cos(theta);-l*sin(theta)];
% O_1to2dot=[0;-thetadot*l*sin(theta);-thetadot*l*cos(theta)];
% O_1to2doubledot=[0;-thetadoubledot*l*sin(theta)-thetadot^2*l*cos(theta);...
%     -thetadoubledot*l*cos(theta)+thetadot^2*l*sin(theta)];
% 
% P_2=[r*cos(phi);r*sin(phi);0];    
% P_2dot=[-r*phidot*sin(phi);r*phidot*cos(phi);0];
% P_2doubledot=[-r*phidoubledot*sin(phi)-r*phidot^2*cos(phi);
%     r*phidoubledot*cos(phi)-r*phidot^2*sin(phi);0];
% omega=[0;0;phidot];

%cross product matrix of rotation rate
% omega_hat=[0 -omega(3) omega(2);
%            omega(3) 0 -omega(1);
%            -omega(2) omega(1) 0];
% 
% theta=theta0+omega(3)*t;
% 
% R_1to2=[cos(theta) -sin(theta) 0;
%       sin(theta) cos(theta) 0;
%        0 0 1];

 
out=[thetadot; thetadoubledot; phidot; phidoubledot];
end