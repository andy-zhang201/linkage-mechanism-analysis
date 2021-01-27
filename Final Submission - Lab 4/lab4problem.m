%%% MIE301 Lab 4
%% Starter Code - Calculate mechanism motion and plot it
close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
xmin= -.04; % leftmost window edge
xmax= .04;  % rightmost window edge
ymin= -.03; % bottom window edge
ymax= .09;  % top window edge

% Link Parameters
theta2_steps = 50; % number of theta2 configuration steps to calculate along mechanism rotation
max_rotation_theta2 = 360 *pi/180; % rotation limit of theta2, radians
theta2 = linspace(0,max_rotation_theta2,theta2_steps); % link 2 rotation into 'increments' number of angles
R = 0.02;          % link #2 length R, cm
L = 0.05;          % link 3 length L, cm
a = 0.015;          % offset, cm
W = 2*9.81;  %weight in N

% set up figure
figure(1);                         %create new figure
% set(1,'WindowStyle','Docked')    %dock the figure

% calculate mechanism motion and plot it
for i=1:theta2_steps                             % step through motion of the mechanism
    hold off;    
    % Draw Link 2:
    Ax(i) = 0;                                   % pivot point of link 2 position
    Ay(i) = 0;                                   % pivot point of link 2 position
    Bx(i) = R*cos( theta2(i) );                  % point B position
    By(i) = R*sin( theta2(i) );                  % point B position
    plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); % draw the link from A to B for the current configuration
    hold on;
    grid on;    
    % Draw Base Pivot:
    recsz = 0.005;                                  % size of drawn base pivot
    plot([0,recsz],[0,-recsz],'r');                 % draw base pivot
    plot([0,-recsz],[0,-recsz],'r');                % draw base pivot
    plot(0,0,'ro','MarkerFaceColor','w');           % draw a small circle at the base pivot point
    plot(Bx(i), By(i), 'bo','MarkerFaceColor','w'); % draw a small circle at B
    text(Bx(i)+.002, By(i), 'B','color','b');       % label point B    
    % Calculate D here:
    theta3(i) = acos( (a-R*cos(theta2(i)))/L );
    D(i) = R*sin(theta2(i)) + L*sin(theta3(i)); 
    % D(i) = R*sin(theta2(i)) + sqrt( L^2 - ( a-R*cos(theta2(i)) )^2 );
    % Draw Link 3:
    plot(a,D(i),'mo','MarkerFaceColor','w');           % draw a small circle at the piston pivot point
    plot( [Bx(i), a], [ By(i), D(i)],'Color','m','LineWidth',3 );   % draw the link from B to C for the current configuration
    % Draw piston:
    rectangle('position',[ a-.0025, D(i)-.005, .005, .01 ] );       % draw the piston itself [x y w h]
    % Figure properties 
    hold on;                            % draw on top of the current figure
    xlabel('x (m)', 'fontsize', 15);    % axis label
    ylabel('y (m)', 'fontsize', 15);    % axil label 
    grid off;                           % add a grid to the figure
    axis equal;                         % make sure the figure is not stretched
    title('Lab4 - Starter');            % add a title to the figure
    axis( [xmin xmax ymin ymax] );      % figure axis limits
    pause(0.009);                       % wait to proceed to next configuration, seconds
    
end


%% Part I - Quasi-Static Force Analysis
%  Part I, a & b
clear all
close all
clc

theta2_steps = 50; 
theta2 = linspace(0, 2*pi, theta2_steps);
R = .02;          % link #2 length R, m
L = .05;          % link 3 length L, m
intervals = 3;
a=[0, 0.01, 0.015];
weight = 2*9.81;
for j=1:intervals
    for i=1:theta2_steps
        theta3(i) = acos( (a(j)-R*cos(theta2(i)))/L );        
%         Formula for M2 in terms of weight
        M2(i) = weight/sin(theta3(i))*R*sin(theta3(i)-theta2(i));
    end

end

[maxMom, maxMomInd] = max(M2);
maxtheta2 = theta2(maxMomInd)


%% Part I, C: Work Calculations 
weight = 2*9.81;

% stroke length
stroke = sqrt((R+L)^2-a(3)^2)-sqrt((R-L)^2-a(3)^2)

% Work Done to lift
W_lift = stroke*weight

MotorWork = 0;  % Initializing accumulator
minTheta = acos(a(3)/(R+L));
maxTheta = acos(a(3)/(R-L));

indMin = find(theta2>minTheta, 1, "first")
indMax = find(theta2<(minTheta + pi), 1, "last")

theta2(indMin);
theta2(indMax);

for i = indMin:(indMax-1)     % Numerical Integral
    MW_partial = (theta2(i+1)-theta2(i))*weight/sin(theta3(i+1))*0.02*sin(theta3(i+1)-theta2(i+1));
    MotorWork = MotorWork + MW_partial;
end
MotorWork


%% Part I, D: Plot torque for several offset values a = 0 to 2cm
clear all
close all
clc

theta2_steps = 50; 
theta2 = linspace(0, 2*pi, theta2_steps);
R = .02;          % link #2 length R, cm
L = .05;          % link 3 length L, cm
weight = 2*9.81;
intervals = 4;
a=linspace(0,0.02,intervals);
% a=[0, 0.01, 0.015];
weight = 2*9.81;

figure(2);

for j=1:intervals
    hold on;
    for i=1:theta2_steps
        theta3(i) = acos( (a(j)-R*cos(theta2(i)))/L );
        D(i) = R*sin(theta2(i)) + L*sin(theta3(i)); 
        M2(i) = weight/sin(theta3(i))*R*sin(theta3(i)-theta2(i));
    end
    plot(theta2, M2);

end
legend("a=0","a=0.007","a=0.013","a=0.02","NumColumns",2)
title('Lab4 - Work done by motor');
xlabel('Theta2 (rad)', 'fontsize', 15);   % axis label
ylabel('M2 (Nm)', 'fontsize', 15);   % axil label


%% Part II,  Dynamic Force Analysis
%  Part II, e & f

clear all
close all
clc

theta2_steps = 100; 
theta2 = linspace(0, 2*pi, theta2_steps);
R = 0.02;          % link #2 length R, cm
L = 0.05;          % link 3 length L, cm
a = 0.015;

rpm_steps = 5;
theta2_dot = linspace(10,100,rpm_steps); % rotation rate, rpm

for j=1:rpm_steps
    t_rev = 60/theta2_dot(j); %Minutes per revolution * 60 seconds
    time = linspace(0, t_rev, theta2_steps);
    %%% as rpm changes, time step changes too! 
    %%% you need timeStep(j) to calculate velocity and acceleration!
    % ∆t = time full revolution/theta2_steps
    timeStep(j) = t_rev/theta2_steps;   
     
    %For j=1, rpm =10
    for i=1:theta2_steps 
        theta3(i) = acos( (a-R*cos(theta2(i)))/L );
        D(i) = R*sin(theta2(i)) + L*sin(theta3(i));
        
        %%% velocity calculations - numerical method
        if i<2
            vel(i,j) = 0; % Velocity is set to 0 when theta2 is 0.
        else
            %When we have more than 1 position, calculate velocity using
            %difference in position over 1 timestep.
            vel(i,j) = (D(i)-D(i-1))/(timeStep(j));  
        end
        %%% acceleration calculations - numerical method
        if i<3
            acc(i,j) = 0; % The change in velocity from the 3rd step onwards is considered. 
        else
            acc(i,j) = (vel(i,j)-vel(i-1,j))/(timeStep(j)); 
        end
        
        %%% Torque calculations
        weight(i, j) = 2*9.81;     
        F_Inertial(i,j) = 2*acc(i,j);

        M2(i,j) = ((weight(i,j)+F_Inertial(i,j))/sin(theta3(i)))*R*sin(theta3(i)-theta2(i));
    end    
end
figure(3);
plot(theta2,M2);
legend('rpm=10','rpm=32.5','rpm=55',"rpm=77.5","rpm=100","NumColumns",2)
xlabel('Theta2 (radians)', 'fontsize', 15); 
ylabel('Moment from motor (Nm)', 'fontsize', 15);   
title('Lab4 - Part II f ');    
vel
acc
M2;

%% Part II, g:

for i=1:5
    [A,B] = max(abs(acc(:,i)));
    maxAccel(i) = A.*sign(acc(B,i));
    
end
maxAccel







%% End

















