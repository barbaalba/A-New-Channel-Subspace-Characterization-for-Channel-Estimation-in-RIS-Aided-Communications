%% To search for the optimum intial points 
%% at the bottom, it includes also an AZ-EL graph
clear; clc; close all;
%% Initialize the parameters
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
RISspacing = 1/4;
M_H = 8;
M_V = 8;
%% Specify the search ranges
azref = -pi/2:pi/180:pi/2;   
elref = -pi/2:pi/180:pi/2;  
elref = deg2rad(0); % one reference angle only
azref = deg2rad(90); % one reference angle only

%% Collects the dimensionality
num = zeros(length(elref),length(azref));
for j = 1:length(elref)
    for i = 1:length(azref)
        [theta,phi,num(j,i)] = UPA_BasisElupnew(M_V,M_H,RISspacing,RISspacing,azref(i),elref(j));  
    end
end

%% Find the optimum combinations
maxval = max(num,[],'all');
[x,y] = find(num == maxval);
elcandidate = elref(x);
elcandidate = rad2deg(elcandidate);
azcandidate = azref(y);
azcandidate = rad2deg(azcandidate);

%% plot dimensionalities for all AZ-EL values
if length(elref)>1 && length(azref)>1
    [X,Y] = meshgrid(azref,elref);
    figure;
    surf(rad2deg(X),rad2deg(Y),num);
    view(-164.3,24.5);
    colorbar;
    xlabel('Azimuth','FontSize',20);
    ylabel('Elevation','FontSize',20);
    title(M_H+ "\times" +M_V);
    % Figure Representation
    ax = gca; % to get the axis handle
    ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
    ax.Position = [0.1 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                               % the figure border
    ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                      % axis
    fig = gcf;
    set(fig,'position',[60 50 1100 800]);

    %% Plot dimensionalities in 2D AZ-EL
    figure;
    imagesc(rad2deg(azref),rad2deg(elref),num);
    xlabel('Azimuth','FontSize',20);
    ylabel('Elevation','FontSize',20);
    title(M_H+ "\times" +M_V);
    colorbar;
end

%% AZ-EL Graph 
figure('DefaultAxesFontSize',20,'defaultLineLineWidth',2,'defaultAxesTickLabelInterpreter','latex')
x = linspace(-90,90,360);
fn = fieldnames(phi);
for i = 1:length(theta)
    y = theta(i);
    plot(x,repelem(rad2deg(y),length(x)),'b'); % blue lines
    hold on;
    for j = 1:length(phi.(fn{i}))
        plot(rad2deg(phi.(fn{i})(j)),rad2deg(y),...
            'r','Marker','x','MarkerSize',15); % red cross
    end
end
xlim([-90,90]);
ylim([-90,90]);
xlabel('Azimuth','Interpreter','latex');
ylabel('Elevation','Interpreter','latex');
xticks(-90:20:90);
xticklabels({'$-90$','$-70$','$-50$','$-30$','$-10$','$10$','$30$','$50$','$70$','$90$'});
yticks(-90:20:90);
yticklabels({'$-90$','$-70$','$-50$','$-30$','$-10$','$10$','$30$','$50$','$70$','$90$'});
grid on;
% Figure Representation
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
