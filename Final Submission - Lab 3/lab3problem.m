%%% MIE301 Lab 3: Optimizing mechanisms
%%
close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
xmin= -20; % leftmost window edge
xmax= 70;  % rightmost window edge
ymin= -25; % bottom window edge
ymax= 75;  % top window edge

%% Link Parameters
increments = 50; % number of theta2 configuration steps to calculate along mechanism rotation %%%% YOU MAY WANT TO CHANGE THIS

max_rotation_theta2 = 360 *pi/180; % rotation limit of theta2, radians
theta2 = linspace(0,max_rotation_theta2,increments); % link 2 rotation into 'increments' number of angles

r1 = 19;          % link 1 length, cm
r2 = 10;          % link 2 length, cm
r3 = 25;          % link 3 length, cm
r4 = 25;          % link 4 length, cm

a = 3.5 ; % path straightness tolerance distance
%% Rotation time
theta2_dot = 60;                        % rotation rate, rpm
t_rev = 60/theta2_dot;                  % rotation time limit, seconds
time = linspace(0,t_rev,increments);    % link 2 rotation into 'increments' number of times

%% set up figure
figure(1);                         %create new figure
% set(1,'WindowStyle','Docked')      %dock the figure

%% 4-bar mechanism geometric constants: See textbook section 4.3.3 for derivation. (Eq. 4.3-54)
h1 = r1/r2;
h2 = r1/r3;
h3 = r1/r4;
h4 =(-r1^2-r2^2-r3^2+r4^2)/(2*r2*r3);
h5 =(r1^2+r2^2-r3^2+r4^2)/(2*r2*r4);

%% define extension vector here; to enable cycling through each extension value
extensionVector = linspace(10,40,5);   % link 3 extension length
                        % define it as a vector to replace "extensionScalar"
                        % and write an outter loop to cycle through each extension value
                        % Hint: after definig the outter loop, Cx(i,1) and Cy(i,1) can be modified!

%% calculate mechanism motion and plot it 
B=size(extensionVector);
for j=1:B(2)   
    for i=1:increments                        % step through motion of the mechanism
        hold off;

        %  geometric calculations (book eq. 4.3-56 to 4.3-62):
        d = -h1 +(1-h3)*cos(theta2(i)) +h5;
        b = -2*sin(theta2(i));
        e = h1 -(1+h3)*cos(theta2(i)) +h5;
        a_a = -h1 +(1+h2)*cos(theta2(i)) +h4;
        c = h1 -(1-h2)*cos(theta2(i)) +h4;

        theta3_1(i) = 2*atan(((-b-(b^2-4*a_a*c)^0.5)/(2*a_a))); %calculate angle of link 3 (eq. 4.3-64)
        theta4_1(i) = 2*atan(((-b-(b^2-4*d*e)^0.5)/(2*d)));          

        % Link Coordinates calculations:   
        Ax(i) = 0;                                      % pivot point of link 2 position
        Ay(i) = 0;                                      % pivot point of link 2 position
        Bx(i) = r2*cos( theta2(i) );                    % point B position
        By(i) = r2*sin( theta2(i) );                    % point B position
        Cx(i,j) = Bx(i) + (r3+extensionVector(j) )*cos( theta3_1(i) );             % point C position. We're going to store these for each extension value (indexed by j)
        Cy(i,j) = By(i) + (r3+extensionVector(j) )*sin( theta3_1(i) );             % point C position
        Dx(i) = r1 + r4*cos( theta4_1(i) );             % point D position
        Dy(i) = r4*sin( theta4_1(i) );                  % point D position        
        
        if (true) % put true here to draw the mechanism each time, false to skip drawing that for increased simulation speed
            plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); % draw link2
            hold on;            
            plot( [Bx(i), Dx(i)], [ By(i), Dy(i)],'Color','b','LineWidth',3 ); % draw link3
            plot( [Dx(i), Cx(i,j)], [ Dy(i), Cy(i,j)],'Color','b','LineWidth',3 ); % draw link3 extension to point C
            plot( [r1, Dx(i)], [ 0, Dy(i)],'Color','m','LineWidth',3 ); % draw link4
            
            % Draw Base Pivots:
            recsz = 2.5;                                        % size of drawn base pivot
            plot([0,recsz],[0,-recsz],'r');                     % draw base pivot for link2
            plot([0,-recsz],[0,-recsz],'r');                    % draw base pivot for link2
            plot(0,0,'ro','MarkerFaceColor','w');               % draw a small circle at the base pivot point
            plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     % draw a small circle at B
            text(Bx(i)+0.9, By(i), 'B','color','b');            % label point B

            plot([r1,r1+recsz],[0,-recsz],'r');                 % draw base pivot for link4
            plot([r1,r1-recsz],[0,-recsz],'r');                 % draw base pivot for link4
            plot(r1,0,'ro','MarkerFaceColor','w');              % draw a small circle at the base pivot point
            plot(Dx(i), Dy(i), 'bo','MarkerFaceColor','w');     % draw a small circle at D
            text(Dx(i)+0.9, Dy(i), 'D','color','b');            % label point D
            text(Cx(i,j)+0.9, Cy(i,j), 'C','color','b');        % label point C
            
            xlabel('x (cm)', 'fontsize', 15);   % axis label
            ylabel('y (cm)', 'fontsize', 15);   % axil label 
            title('Lab3 - starter');            % add a title to the figure
            axis equal;                         % make sure the figure is not stretched
            grid on;
            axis( [xmin xmax ymin ymax] );      % figure axis limits            
        end
        pause(0.01);                            % wait to proceed to next configuration, seconds
    end        
%     pause(.01); %pause after drawing the current and previous paths, seconds



%% plot traveled Path of C on the same figure showing mechanism motion, figure(1)
 %step through all j iterations to plot the straight segments
%  Implemented after I implement j iterations in Cx, Cy
    [lowest_y, index_lowest] = min(Cy(:,j)); % height of lowest point on path
    firstLimit(j) = find(Cy(:,j)<lowest_y+a,1,'first');
    secondLimit(j) = find(Cy(:,j)<lowest_y+a,1,'last');
    Cx_first(j) = Cx(firstLimit(j),j);
    Cx_last(j) = Cx(secondLimit(j),j);
    g(j) = Cx_first(j)-Cx_last(j);

    plot(Cx(:,j),Cy(:,j),"k-")
        
%     pause(0.5);                            % wait to proceed to next configuration, seconds

end
%% part a - Plot your optimization plot, g as a function of the extension length. 
figure(2)
xlabel('x (cm)', 'fontsize', 15);   % axis label
ylabel('y (cm)', 'fontsize', 15);   % axil label 
axis equal;                         % make sure the figure is not stretched
grid on;

axis([xmin xmax ymin ymax] );      % figure axis limits  
plot(extensionVector,g,'-bo');

%% part b - Plot the optimized trace; with the straight line region precisely marked. 
%  Include the optimized mechanism links. 
figure(3)

xlabel('x (cm)', 'fontsize', 15);   % axis label
ylabel('y (cm)', 'fontsize', 15);   % axil label 
title('Lab3 - Optimized Trace');            % add a title to the figure
axis equal;                         % make sure the figure is not stretched
grid on;
axis([xmin xmax ymin ymax] );      % figure axis limits  
hold on;
% Plotting the optimized trace path
[maxG,maxGInd] = max(g);
% Highlight
plot(Cx(firstLimit(maxGInd):secondLimit(maxGInd),maxGInd),Cy(firstLimit(maxGInd):secondLimit(maxGInd),maxGInd),"-r",'LineWidth',3)
plot(Cx(:,maxGInd),Cy(:,maxGInd),"k-")

% Plot trace of mechanism
% for i=1:increments
i=1;
plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); % draw link2
hold on;            
plot( [Bx(i), Dx(i)], [ By(i), Dy(i)],'Color','b','LineWidth',3 ); % draw link3
plot( [Dx(i), Cx(i,maxGInd)], [ Dy(i), Cy(i,maxGInd)],'Color','b','LineWidth',3 ); % draw link3 extension to point C
plot( [r1, Dx(i)], [ 0, Dy(i)],'Color','m','LineWidth',3 ); % draw link4

% Draw Base Pivots:
recsz = 2.5;                                        % size of drawn base pivot
plot([0,recsz],[0,-recsz],'r');                     % draw base pivot for link2
plot([0,-recsz],[0,-recsz],'r');                    % draw base pivot for link2
plot(0,0,'ro','MarkerFaceColor','w');               % draw a small circle at the base pivot point
plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     % draw a small circle at B
text(Bx(i)+0.9, By(i), 'B','color','b');            % label point B

plot([r1,r1+recsz],[0,-recsz],'r');                 % draw base pivot for link4
plot([r1,r1-recsz],[0,-recsz],'r');                 % draw base pivot for link4
plot(r1,0,'ro','MarkerFaceColor','w');              % draw a small circle at the base pivot point
plot(Dx(i), Dy(i), 'bo','MarkerFaceColor','w');     % draw a small circle at D
text(Dx(i)+0.9, Dy(i), 'D','color','b');            % label point D
text(Cx(1,maxGInd)+0.9, Cy(1,maxGInd), 'C','color','b');        % label point C
% end
%% part c - Input link angle?

theta2_dot_r = theta2_dot*2*pi/60;

% Find the time at which first limit and second limits are reached. Equal
% to i.
secondsPerIncrement = t_rev/increments;

time1 = firstLimit(maxGInd)*secondsPerIncrement;
time2 = secondLimit(maxGInd)*secondsPerIncrement;

% Find theta in degrees by multiplying time with angular velocity
thetaLimit1 = theta2_dot_r*time1*180/pi;
thetaLimit2 = theta2_dot_r*time2*180/pi;



%% part d - Average speed of the trace point through the straight region? 

% speed = total distance/total time

tot_time = time2-time1;
tot_dist = g(maxGInd);

speed = tot_dist/tot_time;


%% part e - How does the tolerance (a) affect your optimized path length? 

% Increasing a will increase the optimized path length.




