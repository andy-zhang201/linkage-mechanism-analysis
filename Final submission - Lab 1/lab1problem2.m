%%% MIE301 Lab 1
%%
close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
% using axis()
xmin= -8;  % leftmost window edge
xmax= 8;   % rightmost window edge
ymin= -15; % bottom window edge
ymax= 15;  % top window edge

%% Define the link and motion Parameters
stepsize = 30 *pi/180;                          % number of configuration steps to calculate along mechanism rotation
max_rotation_theta2 = 300 *pi/180;              % rotation limit of theta2, radians
theta2 = 0 : stepsize : max_rotation_theta2;    % link 2 rotation into 'increments' number of angles
phi = 33 * pi/180;          % Fixed angle in link 2 in radians. This should be the last 2 digits of my student number.
length2 = 2;                                    % link #2 length
myVec1 = [];
myVec2 = [];


%% set up figure
figure(1);                         %create new figure
set(1,'WindowStyle','Docked')      %dock the figure

%% calculate mechanism motion and plot it
for i=1:length(theta2)                            % step through motion of the mechanism
       
    %%I comment out the hold off command in order to plot the instances of
    %%motion from each mechanism. I wrote a code file for each mechanism that I had to animate.
    %%I made the figures for problem two by just commenting out the hold off command.
    %%This code animates the motion for mechanism 2. 
    
    %%
    
    %hold off;                                     % erase the previously plotted data each cycle
    Ax(i) = 0;                                    % pivot point of link 2 position
    Ay(i) = 0;                                    % pivot point of link 2 position
    Bx(i) = length2*cos( theta2(i) );             % point B position
    By(i) = length2*sin( theta2(i) );             % point B position
    plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); %draw the link from A to B for the current configuration
    
    %Save link AB graph
    hold on;
    
    
    
    basePivotSize = .5;                                 %size of drawn base pivot
   
    plot([0,basePivotSize],[0,-basePivotSize],'r');     %draw base pivot
    plot([0,-basePivotSize],[0,-basePivotSize],'r');    %draw base pivot
    plot(0,0,'ro','MarkerFaceColor','w');               %draw a small circle at the base pivot point
    plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     %draw a small circle at B    
    
    hold on;                            % draw on top of the current figure
    xlabel('x (cm)', 'fontsize', 15);   % axis label
    ylabel('y (cm)', 'fontsize', 15);   % axil label 
    grid on;                            % add a grid to the figure
    axis equal;                         % make sure the figure is not stretched
    title('Lab1 - Problem 1');          % add a title to the figure
    axis( [xmin xmax ymin ymax] );      % set figure axis limits
    
    %% add new code here:   
    
    %%Configuration Plot Mek A
    %Define variables
    magnitude_BC = 3;
    magnitude_BA = 2;
    phi = 33*pi/180;
    
    %Finding components of vector BC using coordinates for vector AB
    BCx = (magnitude_BC/magnitude_BA)*(-cos(phi)*(Bx(i))+sin(phi)*(By(i)));
    BCy = (magnitude_BC/magnitude_BA)*(-sin(phi)*(Bx(i))-cos(phi)*By(i));
    
    %Finding C coordinates by adding vector BC to vector AB
    Cx = Bx(i)+BCx;
    Cy = By(i)+BCy;
    
    myVec1 = [myVec1,Cx];
    myVec2 = [myVec2,Cy];
    
    %Plotting Link BC
    plot([Bx(i),Cx],[By(i),Cy],"b-","LineWidth",3);
    
    %Points
    plot(Cx,Cy,"bo","MarkerFaceColor","w");             %small circle C
    %text(Cx +0.4, Cy, 'C','color','b');                      %label point C
    %Save link graph
    hold on;   
    
    
    
    %%
    pause(0.1);                       % wait to proceed to next configuration, in seconds
    
end
%labels
text(Cx +0.4, Cy, 'C','color','b');
text(Bx(i)+0.4, By(i), 'B','color','b'); 

plot(Bx(:),By(:),'--g','linewidth',1);  % draw a trace of the path of point B
plot(myVec1(:),myVec2(:),'--m','linewidth',1);

%%





