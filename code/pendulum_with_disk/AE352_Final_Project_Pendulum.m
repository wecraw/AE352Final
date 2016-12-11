clear;clc;close all;
l=20;%pendulum length in meters
r=12;%disk radius in meters
g=9.8;%m/s^2
m=2000;%mass in kilograms of disk
J=m*l^2;%moment of inertia
F_tire=1000;%Newtons
L_tire=.5;%contact area of tire in meters
T_disk=1000;%torque on disk
J_disk=.5*m*r^2;%moment of inertia of disk in Newton-meters
oscillations=15;
Numpeople=30;

t0=0;
theta0=.03;%radians
thetadot0=0.3;%radians/second (give a kick to start the simulation
phi0vec=linspace(0,2*pi,Numpeople);%starting position of each rider
phidot0=1;%radians/second
tlength=50;%number of points for each numerical integration
for w=1:length(phi0vec)%run through each 
    theta0=.03;%radians
    thetadot0=0.3;%radians/second (give a kick to start the simulation
    phi0vec=linspace(0,2*pi,Numpeople)
    phidot0=1;%radians/second
    phi0=phi0vec(w);%current rider to work with
    t_out=[];
    y_out=[];
    %output vectors (appended on each run to show total run of pendulum)
    for i=1:oscillations
        if i==1
        thetamax=acos(1-thetadot0^2*l/(2*g)); %find max theta value (see notes for derivation)
        T=2*pi*sqrt(abs(l/g))*(1+1/16*thetamax^2+11/3072*thetamax^4);%period

        infovec=[theta0 thetadot0 phi0 phidot0];
        %tvec=t0:tstep:(T/2);%this for sure works
        tvec=linspace(t0,T/2,tlength);
        [t, y]= ode45(@thetafunc, tvec,infovec);%go from t=0 to when pendulum completes...
        %half an oscillation (back at the bottom)
        t_out=t;
        y_out=y;
        tf=T/2;
        else
        %load in new information from last ODE 45 run
        thetadot0=y(size(y,1),2);
        phi0=y(size(y,1),3);
        phidot0=y(size(y,1),4);

        v=thetadot0*l;%velocity at bottom
        t_contact=L_tire/v; %approximate time of contact between ride and tire
        torque_tire=F_tire*l;%torque exerted by the tire
        deltaH=torque_tire*t_contact;%change in angular momentum
        thetadotf=thetadot0+deltaH/J;%angular velocity after tire

        thetamax=acos(1-thetadotf^2*l/(2*g));
        T=2*pi*sqrt(abs(l/g))*(1+1/16*thetamax^2+11/3072*thetamax^4);

        infovec=[theta0 thetadotf phi0 phidot0];
        %tvec=tf:tstep:(tf+T/2);%this for sure works
        tvec=linspace(tf,T/2+tf,tlength);
        [t, y]= ode45(@thetafunc, tvec,infovec);%go from end of last run...
        %to when pendulum completes half an oscillation (back at the bottom)
        tf=tf+T/2;
        %append outputs to outputs from previous runs to plot later
        t_out=[t_out;t];
        y_out=[y_out;y];
        end
    end


    %outputs
    thetaout=y_out(:,1).';
    thetadotout=y_out(:,2).';
    phiout=y_out(:,3).';
    phidotout=y_out(:,4).';
    %rotation matrix from disk to base frame
    R0to1=[ones(1,length(thetaout)) zeros(1,length(thetaout)) zeros(1,length(thetaout));...
        zeros(1,length(thetaout)) cos(thetaout) -sin(thetaout);...
        zeros(1,length(thetaout)) sin(thetaout) cos(thetaout)];
    P1=[r*cos(phiout); r*sin(phiout); zeros(1,length(thetaout))];%point in disk frame
    %initialize storage vectors
    P0=[];
    P1=[];
    P2=[];
    Origin=[];
    for i=1:length(thetaout)%runs through all theta values for pendulum
        %origin transformation to base frame
        O_0to1=[0;l*sin(thetaout(i));-l*cos(thetaout(i))];
        %rotation matrix to base frame
        R0to1=[1 0 0;...
        0 cos(thetaout(i)) -sin(thetaout(i));...
        0 sin(thetaout(i)) cos(thetaout(i))];
        P1=[r*cos(phiout(i)); r*sin(phiout(i)); 0];
        Origin=[Origin,O_0to1];%store different frame 1 origin values (used for plotting pendulum motion)
        P0=[P0,O_0to1+R0to1*P1];%store information about point in base frame
        %^above outputs x,y,z in rows 1, 2, 3
    end
    %initialize and store point vectors for each person
    if w==1%if first time, create storage arrays
        P0array=P0;
        OriginArray=Origin;
    else%for each rider, store position and relative origin information
        P0array(:,:,w)=P0;
        OriginArray(:,:,w)=Origin;
    end
end
%plot theta output
figure;
plot(t_out,y_out(:,1));
xlabel('Time, t, seconds')
ylabel('Angle, \theta, radians')
title('\Theta vs. t, Paul DeTrempe, AE 352 Pirate Ship Model')

%plot 3-Space
figure;
xlabel('x');
ylabel('y');
zlabel('z');

    for j=1:size(Origin,2)
        personvector=[];%used for storing points of people
        %axis limits to make it look pretty
        xlim([-25 25]);
        ylim([-25 25]);
        zlim([-25 25]);
        grid on;
        
        for count=1:size(P0array,3)%run through each rider
        %plot point on disk for each rider
        person=plot3(P0array(1,j,count),P0array(2,j,count),P0array(3,j,count),'ro');
        personvector=[personvector,person];%store plotted points for later deletion
        end
        hold on;
        %plot a line from (0,0,0) to the frame 1 origin point. this is the
        %pendulum
        pendulum=plot3([0, Origin(1,j)],[0,Origin(2,j)],[0,Origin(3,j)],'g');
        pause(.01);
        hold on;
        %delete previous plotted points and pendulum
        delete(pendulum);  
        delete(personvector);
    end
% end
