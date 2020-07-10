% Simple finite difference solution to the 2D wave equation. Gradients in
% time and space are both calculated using a 2nd-order accurate central
% difference scheme. Adapted from 1D solution provided by Dr Bradley
% Treeby (Ultrasound in Medicine Course - MPHYx900 - UCL).

% clear all variables in the workspace, close all active figures and clear
% command window

clc
close all
clear all

% ----------
% Question 3
% ----------

% set the literals (hard-coded numbers used in the script) and initial
% conditions

% number of grid points in x and y directions
Nx = 100;
Ny = 100;

% grid spacing (m) in x and y directions
dx = 1e-3;  
dy = 1e-3;

% sound speed (m/s)
c0 = 1500;

% number of time steps
Nt = 100;

% create the grid axis
x = (1:Nx)*dx; 
y = (1:Ny)*dy;

% set the position of the source
x_pos = (Nx/2)*dx;
y_pos = (Ny/2)*dy; 

% set the initial pressure to be a Gaussian
variance = (2*dx)^2;
gaussian_x = exp( -(x - x_pos).^2 / (2 * variance) );
gaussian_y = exp( -(y - y_pos).^2 / (2 * variance) );
p_n = gaussian_x' * gaussian_y; 

% set the size of the time step to its maximum within the stability limit, 
% where in 2D the stability limit is given by dt <= dx/(sqrt(2)*c0)
dt = dx/(sqrt(2)*c0);

% set pressure at (n - 1) to be equal to pressure at n (this implicitly sets the
% initial particle velocity to be zero)
p_nm1 = p_n;

% preallocate the pressure at (n + 1) (this is updated during the time
% loop)
p_np1 = zeros(size(p_n));

% open a new figure in a maximised window (to facilitate saving images)
figure('units','normalized','outerposition',[0 0 1 1])

%to initialise vectors containing max and min z values for each time step
%(to facilitate setting appropriate z axis correctly later)
min_vector=zeros(1,Nt);
max_vector=zeros(1,Nt);

%to initialise vector containing value of pressure field at grid position
% (40,40), which is used in Question 6
p_out=zeros(1,Nt);

% To set line widths (for 2D plots) and font size (for all plots)
LW=1.5;
fs=20;

% calculate pressure in a loop through the time steps
for n = 1:Nt
    
    % -----------
    % CALCULATION
    % -----------
    
    % calculate the new value for the pressure, leaving edge values as
    % zero
    
    %let T = time derivative terms of 2D solution
    T = 2*p_n(2:end-1,2:end-1) - p_nm1(2:end-1,2:end-1);
    %let X = x derivative terms of 2D solution (noting Matlab has 
    %(row,column) --> (Y,X) notation when addressing a matrix location)
    X = p_n(2:end-1,1:end-2) - 2*p_n(2:end-1,2:end-1) + p_n(2:end-1,3:end);
    %let Y = y derivative terms of 2D solution
    Y = p_n(1:end-2,2:end-1) - 2*p_n(2:end-1,2:end-1) + p_n(3:end,2:end-1);
    
    %calculate pressure at n+1 in 2D using the above three terms
    p_np1(2:end-1,2:end-1) = T + ( c0^2 * dt^2 * (X/dx^2 + Y/dy^2) ); 
     
    %copy the value of p at n to p at (n - 1)
    p_nm1 = p_n;
    
    % copy the pressure at (n + 1) to the pressure at n
    p_n = p_np1;
    
    % --------
    % PLOTTING
    % --------
    
    % plot the pressure field (in units of cm and pascals)
    s=surf(x*100,y*100,p_np1);
    
    % set the limits on the z-axis (derived from max_vector and min_vector 
    % at end of loop)
    set(gca, 'ZLim', [-0.3, 0.8]);
    
    % add a title, subtitle, label axes, and colorbar
    title({'2D Solution to the Wave Equation'; ...
        'Time Step Set to Stability Limit';['n = ' ...
        ,num2str(n),'/',num2str(Nt)]},'FontSize',fs+1,'FontWeight','bold')
    xlabel('x / cm','FontWeight','bold','FontSize',fs-1)
    ylabel('y / cm','FontWeight','bold','FontSize',fs-1)
    zlabel('Pressure Field / Pa','FontWeight','bold','FontSize',fs-1)
    colorbar
    
    % formatting gridlines and axes label fontsizes
    grid minor
    ax = gca;
    ax.FontSize = fs-3;
   
    % force the plot to update
    drawnow;
    
    % briefly pause before continuing the loop
    pause(0.1)
    
    % to update max and min pressure values for this time step to overall
    % max_vector and min_vector
    max_vector(n)=max(max(p_np1));
    min_vector(n)=min(min(p_np1));
    
    %to save plots at time steps 10, 50, 100
    if n == 10
        saveas(figure(1),'Plot_Time_Step_10','jpg');
    elseif n==50
        saveas(figure(1),'Plot_Time_Step_50','jpg');
    elseif n ==100
        saveas(figure(1),'Plot_Time_Step_100','jpg');
    end 
    
    %to save the value of the pressure field at (40,40) in vector p_out
    p_out(n)=p_np1(40,40);
    
end

%to compute the overall max and min z values for all time steps, in order 
%to allow for suitable limits to be implemented on z-axis
zmax=max(max_vector); % =0.7650 --> round to 0.8
zmin=min(min_vector); % =-0.2940 --> round to -0.3

% ----------
% Question 5
% ----------

%Increase the time step to 2% above the stability limit
dt = 1.02 * (dx/(sqrt(2)*c0));

%Reset pressure matrix values and open a new figure
p_n = gaussian_x' * gaussian_y; 
p_nm1 = p_n;
p_np1 = zeros(size(p_n));
figure('units','normalized','outerposition',[0 0 1 1])

%Repeat above calculations and plotting from Question 3 with increased time
%step value

for n = 1:Nt
    
    % -----------
    % CALCULATION
    % -----------

    %forming terms for the expressions in the 2D solution
    T = 2*p_n(2:end-1,2:end-1) - p_nm1(2:end-1,2:end-1);
    X = p_n(2:end-1,1:end-2) - 2*p_n(2:end-1,2:end-1) + p_n(2:end-1,3:end);
    Y = p_n(1:end-2,2:end-1) - 2*p_n(2:end-1,2:end-1) + p_n(3:end,2:end-1);
    
    %calculate pressure at n+1 in 2D using the above three terms
    p_np1(2:end-1,2:end-1) = T + ( c0^2 * dt^2 * (X/dx^2 + Y/dy^2) ); 
     
    %update the pressure matrix values
    p_nm1 = p_n;
    p_n = p_np1;
    
    % -----------
    % PLOTTING
    % -----------
    
    % plot the pressure field (in units of cm and pascals)
    surf(x*100,y*100,p_np1);
    
    % set the limits on the z-axis to the same as Question 3
    set(gca, 'ZLim', [-0.3, 0.8]);
    
    % add a title, subtitle, label axes, and colorbar
    title({'2D Solution to the Wave Equation'; ...
        'Time Step Set to Beyond Stability Limit';['n = ' ...
        ,num2str(n),'/',num2str(Nt)]},'FontSize',fs+1,'FontWeight','bold')
    xlabel('x / cm','FontWeight','bold','FontSize',fs-1)
    ylabel('y / cm','FontWeight','bold','FontSize',fs-1)
    zlabel('Pressure Field / Pa','FontWeight','bold','FontSize',fs-1)
    colorbar
    
    % formatting gridlines and axes label fontsizes
    grid minor
    ax = gca;
    ax.FontSize = fs-3;
   
    % force the plot to update
    drawnow;
    
    % briefly pause before continuing the loop
    pause(0.1)
    
    %to save plot when instability begins (at roughly time step 95 by 
    %inspection)
    if n == 95
        saveas(figure(2),'Plot_Time_Step_95_Unstable','jpg');
    end 
    
end

% ----------
% Question 6
% ----------

%To plot vector p_out against time (p_out = pressure at grid position
%(40,40))
time = dt*(1:Nt); %first form a time vector in seconds
figure('units','normalized','outerposition',[0 0 1 1])
plot(time,p_out,'color','b','linewidth',LW)
title({'Pressure at Grid Position (40,40)'; ...
        'Gaussian Initial Pressure Variance = (2dx)^{2}'}, ...
        'FontSize',fs+1,'FontWeight','bold')
xlabel('Time / s','FontWeight','bold','FontSize',fs-1)
ylabel('Pressure / Pa','FontWeight','bold','FontSize',fs-1)

%plot formatting
grid minor
ax = gca;
ax.FontSize = fs-3;
xlim([min(time) max(time)]) %tight x-axis fit
line(xlim,[0 0],'color','k');  % plot x-axis

%to save plot
saveas(figure(3),'Pressure_Field_40_40_Original_Variance','jpg');

% ----------
% Question 7
% ----------

%Repeat analysis with reduced variance of Gaussian initial pressure
variance = dx^2;
gaussian_x = exp( -(x - x_pos).^2 / (2 * variance) );
gaussian_y = exp( -(y - y_pos).^2 / (2 * variance) );
%Reset time step value to stability limit
dt = (dx/(sqrt(2)*c0));
%reset pressure matrix values
p_n = gaussian_x' * gaussian_y; 
p_nm1 = p_n;
p_np1 = zeros(size(p_n));

%to initialise a vector containing value of pressure field at grid position
% (40,40) for reduced variance
p_out_reduced_var=zeros(1,Nt);

%Rerun the simulation
for n = 1:Nt
    
    % -----------
    % CALCULATION
    % -----------

    %forming terms for the expressions in the 2D solution
    T = 2*p_n(2:end-1,2:end-1) - p_nm1(2:end-1,2:end-1);
    X = p_n(2:end-1,1:end-2) - 2*p_n(2:end-1,2:end-1) + p_n(2:end-1,3:end);
    Y = p_n(1:end-2,2:end-1) - 2*p_n(2:end-1,2:end-1) + p_n(3:end,2:end-1);
    
    %calculate pressure at n+1 in 2D using the above three terms
    p_np1(2:end-1,2:end-1) = T + ( c0^2 * dt^2 * (X/dx^2 + Y/dy^2) ); 
     
    %update the pressure matrix values
    p_nm1 = p_n;
    p_n = p_np1;
    
    %to save the value of the pressure field at in vector p_out_reduced_var
    p_out_reduced_var(n)=p_np1(40,40);
    
end

%To plot p_out_reduced_var against time
figure('units','normalized','outerposition',[0 0 1 1])
plot(time,p_out_reduced_var,'color','b','linewidth',LW)
title({'Pressure at Grid Position (40,40)'; ...
        'Gaussian Initial Pressure Variance = (dx)^{2}'},'FontSize',fs+1 ...
        ,'FontWeight','bold')
xlabel('Time / s','FontWeight','bold','FontSize',fs-1)
ylabel('Pressure / Pa','FontWeight','bold','FontSize',fs-1)
%plot formatting
grid minor
ax = gca;
ax.FontSize = fs-3;
xlim([min(time) max(time)]) %tight x-axis fit
line(xlim,[0 0],'color','k');  % plot x-axis
%to save plot
saveas(figure(4),'Pressure_Field_40_40_Reduced_Variance','jpg');

%close all open figures
close all
