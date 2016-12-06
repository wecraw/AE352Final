close all;clear;clc;
l=30;%pendulum length in meters
g=9.8;%m/s^2
m=2000;%mass in kilograms of disk
J=m*l^2;%moment of inertia
F_tire=100000;%Newtons
L_tire=.5;%contact area of tire in meters
oscillations=15;

t0=0;

theta0=0;%radians
thetadot0=.05;%radians/second (give a kick to start the simulation

%output vectors (appended on each run to show total run of pendulum)
t_out=[];
y_out=[];
tstep=0.05;
for i=1:oscillations
    if i==1
    thetamax=acos(1-thetadot0^2*l/(2*g)); %find max theta value (see notes for derivation)
    T=2*pi*sqrt(l/g)*(1+1/16*thetamax^2+11/3072*thetamax^4);%period

    infovec=[theta0 thetadot0];
    [t, y]= ode45(@thetafunc, t0:tstep:(T/2),infovec);%go from t=0 to when pendulum completes...
    %half an oscillation (back at the bottom)
    t_out=t;
    y_out=y;
    tf=T/2;
    else
    thetadot0=y(size(y,1),2);%load in new angular velocity at point of contact
    v=thetadot0*l;%velocity at bottom
    t_contact=L_tire/v; %approximate time of contact between ride and tire
    torque_tire=F_tire*l;%torque exerted by the tire
    deltaH=torque_tire*t_contact;%change in angular momentum
    thetadotf=thetadot0+deltaH/J;%angular velocity after tire
    
    thetamax=acos(1-thetadotf^2*l/(2*g));
    T=2*pi*sqrt(l/g)*(1+1/16*thetamax^2+11/3072*thetamax^4);

    infovec=[theta0 thetadotf];
    [t, y]= ode45(@thetafunc, tf:tstep:(tf+T/2),infovec);%go from end of last run...
    %to when pendulum completes half an oscillation (back at the bottom)
    tf=tf+T/2;
    t_out=[t_out;t];
    y_out=[y_out;y];
    end
end
    


%plot outputs
figure(1);
plot(t_out,y_out(:,2));
xlabel('Time, t, seconds')
ylabel('Anglular velocity, \omega, radians/second')
title('\omega vs. t, Paul DeTrempe, AE 352 Pirate Ship Model')
grid on;
figure(2);
plot(t_out,y_out(:,1));
xlabel('Time, t, seconds')
ylabel('Angle, \theta, radians')
title('\theta vs. t, Paul DeTrempe, AE 352 Pirate Ship Model')
grid on;

%Animation

figure(3);
O = [0 0]; %set origin
axis(gca,'equal'); %set aspect ratio of plot
axis([-40 40 -40 20]); %plot limits
grid on;

w_rect=6;
l_rect=3;

support_length=17;
support_height=35;

for i=1:length(t_out);
    P= l*[sin(y_out(i,1)) -cos(y_out(i,1))]; %mass point
    
    O_circ = viscircles(O,1); %circle to represent joint about which 
                                 %the pendulum oscilates
    
    pend = line([O(1) P(1)],[O(2) P(2)],'LineWidth',2);
    
    support1 = line([O(1) -support_length],[O(2) -support_height],'LineWidth',4);
    
    support2 = line([O(1) support_length],[O(2) -support_height],'LineWidth',4);

    
    rect.vertices=[-w_rect/2 -l_rect/2;w_rect/2 -l_rect/2; w_rect/2 l_rect/2; -w_rect/2 l_rect/2]; 

    rect.faces=[1 2 3 4]; %connect vertices
    
    %rotation matrix for ani
    theta =y_out(i,1);
    R=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    ship = patch(rect,'Vertices',rect.vertices*R+repmat([P(1) P(2)],4,1),'FaceColor',[1 0 0]);
    


     %pendulum
    
    pause(0.001); %updates the plot
    
    
    %refresh display
    if i<length(t_out)
        delete(pend);
        delete(ship);
        delete(O_circ);
        delete(support1);
        delete(support2);
    end
    
    
end




