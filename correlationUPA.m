% This code is written to demonstrate the descending eigenvalues of the 
% spatial correlation of a planar array and mark the estimated dimension
% through our proposed method
close all; clear; clc;
%% Initializaation
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 32; % horizental antennas
M_V = 32; % vertical antennas
M = M_V * M_H; % total number of antennas in URA
d_H = [1/4,1/8]; % Antenna spacing
NUE = 100000; % number of instances
%% Find the DOF through our algorithm
azref = pi/2; elref = 0;
DOF = zeros(length(d_H),1);
for i = 1:length(d_H)
    [~,~,DOF(i)] = UPA_BasisElup(M_V,M_H,d_H(i),d_H(i),azref,elref);
end

%% Spatial Correlation 
RIS_AOA = unifrnd(-pi/3,pi/3,1,NUE);
RIS_EOA = unifrnd(-pi/2,pi/2,1,NUE);
for i = 1:length(d_H)
    R = UPA_Evaluate(lambda,M_V,M_H,RIS_AOA,RIS_EOA,d_H(i),d_H(i));
    R = R * R' / NUE; % Correlation matrix
    R = R + R' / 2; % to fix numerical problem of complex eigenvalues
    v = flip(eig(R));
    v = v(v > 0);
    v = pow2db(v);
    plot(v,LineWidth=2); % plot the eigenvalues in descending order
    hold on;
    plot(DOF(i),v(DOF(i)),'x',MarkerSize=15,LineWidth=2);
end
set(groot,'DefaultAxesFontSize',20,'defaultLineLineWidth',2); 
xlim([0,400]);
grid on;
plot(1:400,repelem(0,400,1),'--k','LineWidth',2);
legend('\lambda/4','\eta, \lambda/4','\lambda/8','\eta, \lambda/8','Uncorrelated Rayleigh');
xlabel('Order of Eigenvalues');
ylabel('Eigenvalues [dB]');
%% Figure Configuration
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
