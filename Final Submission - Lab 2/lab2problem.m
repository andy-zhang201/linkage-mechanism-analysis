%%% MIE301 Lab 2
%%
close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
xmin= -30; % leftmost window edge
xmax= 30;  % rightmost window edge
ymin= -20; % bottom window edge
ymax= 50;  % top window edge

%% Link Parameters
increments = 30; % number of theta2 configuration steps to calculate along mechanism rotation
max_rotation_theta2 = 360 *pi/180; % rotation limit of theta2, radians
theta2 = linspace(0,max_rotation_theta2,increments); % link 2 rotation into 'increments' number of angles. Therefore, angle step size is 12 deg
R = 10;         % link #2 length R, cm
L = 30;         % link 3 length L, cm
a = 5;          % offset, cm

%% set up figure
figure(1);                         %create new figure
set(1,'WindowStyle','Docked')      %dock the figure


%% calculate mechanism motion and plot it

% %Preallocation for speed (optional)
Ax = zeros(1,increments);                         
Ay = zeros(1,increments);                        
Bx = zeros(1,increments);                        
By = zeros(1,increments);                         

for i=1:increments                               % step through motion of the mechanism
    
    % Draw Link 2:
    Ax(i) = 0;                                   % pivot point of link 2 position
    Ay(i) = 0;                                   % pivot point of link 2 position
    Bx(i) = -R*sin( theta2(i) );                  % point B position
    By(i) = R*cos( theta2(i) );                  % point B position
    plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); % draw the link from A to B for the current configuration
    hold on;
    grid on;
    
    % Draw Base Pivot:
    recsz = 1.5;                                        % size of drawn base pivot
    plot([0,1.5],[0,-1.5],'r');                         % draw base pivot
    plot([0,-1.5],[0,-1.5],'r');                        % draw base pivot
    plot(0,0,'ro','MarkerFaceColor','w');               % draw a small circle at the base pivot point
    plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     % draw a small circle at B
    
    % Calculate D here:
    Dx(i)= Ax(i)-a;
    D_length = 30 
%     =sqrt((Bx-offset)^2+height^2)
%   height = sqrt(30^2-(Bx-(offsetfromOrigin))^2). Height is relative to point B.

    Height_DB = sqrt(30^2-(abs(Bx(i)+a))^2);
    Dy(i) = By(i)+Height_DB;
    
    plot(Dx(i),Dy(i),'go');
    plot([Bx(i) Dx(i)], [By(i) Dy(i)],'Color','b','LineWidth',3 );
    
    % Draw Link 3:
   

    % Draw piston:
    width = 10;
    height = 20;
    rectangle("Position",[Dx(i)-width/2,Dy(i)-height/2,width,height]);
    
    hold on;                           % draw on top of the current figure
    %hold off;                          % redraw the current figure
    xlabel('x (cm)', 'fontsize', 15);  % axis label
    ylabel('y (cm)', 'fontsize', 15);  % axil label 
    grid off;                          % add a grid to the figure
    axis equal;                        % make sure the figure is not stretched
    title('Lab 2 Figure 1');          % add a title to the figure
    axis( [xmin xmax ymin ymax] );     % figure axis limits
    %axis equal;                      % make sure the figure is not stretched
    pause(0.01);                       % wait to proceed to next configuration, seconds
end
 
%% Rotation time
%theta2_dot = ; % rotation rate, Hz
%t_rev = ; % rotation time limit, seconds
%time = ; % link 2 rotation into 'increments' number of times

%% do final plotting here:
figure(2);
set(2,'WindowStyle','Docked')      %dock the figure
plot(theta2*180/(2*pi),Dy,"ro");
xlabel('Rotation Angle (degree)', 'fontsize', 15);  % axis label
ylabel('D (cm)', 'fontsize', 15);  % axil label 
grid off;                          % add a grid to the figure
title('Lab 2 Plot 1');          % add a title to the figure
hold on;

figure(3);
set(3,'WindowStyle','Docked')      %dock the figure
time = linspace(0,1,30);
plot(time,Dy,"bo");
xlabel('Time (sec)', 'fontsize', 15);  % axis label
ylabel('D (cm)', 'fontsize', 15);  % axil label 
grid off;                          % add a grid to the figure
title('Lab 2 Plot 2');          % add a title to the figure

% calculate stroke ere:
%Angles for limit positions
[maxD,i] = max(Dy)
maxTheta2 = theta2(i)

[minD,j] = min(Dy)
minTheta2 = theta2(j)

%Time for limit positions
maxTime = time(i)
minTime = time(j)

%Stroke legnth
strokeLength = maxD - minD

%Velocity vector

Vel_D = abs(diff(Dy))
maxVelD = max(Vel_D)

hold on;
% figure(2);
% plot(maxTheta2,maxD,"mo");
% plot(minTheta2,minD,"mo");

