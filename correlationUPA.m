% This code is written to demonstrate the descending eigenvalues of the 
% spatial correlation of a planar array
close all; clear; clc;
%% Initializaation
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 32; % horizental antennas
M_V = 32; % vertical antennas
M = M_V * M_H; % total number of antennas in URA
RISspacing = [1/4,1/8]; % Antenna spacing
DOF = [197,49]; % For 32*32 
NUE = 100000; % number of instances
%% Spatial Correlation 
RIS_AOA = unifrnd(-pi/2,pi/2,1,NUE);
RIS_EOA = unifrnd(-pi/2,pi/2,1,NUE);
for i = 1:length(RISspacing)

    R = UPA_Evaluate(lambda,M_V,M_H,RIS_AOA,RIS_EOA,RISspacing(i),RISspacing(i));
    R = R * R' / NUE; % Correlation matrix
    R = R + R' / 2;
    v = flip(eig(R));
    v = v(v > 0);
    v = pow2db(v);
    plot(v,LineWidth=2); % plot the eigenvalues in descending order
    hold on;
    plot(DOF(i),v(DOF(i)),'x',MarkerSize=15,LineWidth=2);
end
xlim([0,400]);
grid on;
plot(1:400,repelem(0,400,1),'--k','LineWidth',2);
legend('\lambda/4','\eta, \lambda/4','\lambda/8','\eta, \lambda/8','Uncorrelated Rayleigh');