%% This code is written to disprove the DFT codebook matrix
close all; clear; clc;
%% Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
%Number of antennas in verical and horizental
M = 4;
%Antenna spacing in wavelengths
interAntennaSpacing = 1/4;
%Prepare to plot colors on a sphere
N = 500;
[X,Y,Z] = sphere(N);

%% Prepare DFT Codebook
W = DFTBookBuild(M,M);

%% For each precoder plot the sphere beamforming plot
for i = 1:M^2
%Prepare to compute channel gains on the sphere
gainMap = zeros(size(X));

    %Go through all azimuth and elevation angles
    for n = 1:size(X,1)
        for m = 1:size(X,2)
        
            %Compute received power according to (7.28) in "Massive MIMO networks"
            [phi2,theta2] = cart2sph(X(n,m),Y(n,m),Z(n,m));
            g = abs(UPA_Evaluate(lambda,M,M,phi2,theta2,interAntennaSpacing,interAntennaSpacing)'*W(:,i)).^2;
            if g > 1
                gainMap(n,m) = g;
            end
            if X(n,m) < 0
                gainMap(n,m) = 0;
            end
       
        end
    end
    % Plot beamforming gain
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 

    figure;
    surf(X,Y,Z,pow2db(gainMap),'EdgeColor','none');
    xlabel('$x$','Interpreter','Latex');
    ylabel('$y$','Interpreter','Latex');
    zlabel('$z$','Interpreter','Latex');
    caxis([-20 pow2db(M^2)+5]); % limit the color map
    colormap(flipud(hot));
    hBar = colorbar;
    set(hBar, 'TickLabelInterpreter', 'latex');
    axis equal;
    set(gca,'color',[0.9 0.9 0.9]);
    set(gca,'fontsize',18);
    view(122,30);
    %xlim([0 1]);
    hold on;
    % Plot the dot lines around the sphere
    varphiAngles = linspace(-pi/2,pi/2,100);
    x_circ = cos(varphiAngles);
    y_circ = sin(varphiAngles);
    plot3(x_circ,y_circ,zeros(size(x_circ)),'k:','LineWidth',2);
end

%Compute interference gains on the sphere
gainMap = zeros(size(X));
[theta,phi] = UPA_BasisEl(lambda,M,M,interAntennaSpacing,interAntennaSpacing);
W = UPA_Codebook(lambda,theta,phi,M,M,interAntennaSpacing,interAntennaSpacing);
for i = 1:size(W,2)
    for n = 1:size(X,1)
        for m = 1:size(X,2)
        
        %Compute received power according to (7.28) in "Massive MIMO networks"
            [phi2,theta2] = cart2sph(X(n,m),Y(n,m),Z(n,m));
            gainMap(n,m) = M^2 * abs(UPA_Evaluate(lambda,M,M,phi2,theta2,interAntennaSpacing,interAntennaSpacing)'*W(:,i)).^2;
           
        end
    end

figure;
surf(X,Y,Z,10*log10(gainMap),'EdgeColor','none');
xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
caxis([-40 pow2db(M^2)]);
colormap(flipud(hot));
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');
axis equal;
set(gca,'color',[0.9 0.9 0.9]);
set(gca,'fontsize',18);
view(122,30);
%xlim([0 1]);
hold on;
plot3(x_circ,y_circ,zeros(size(x_circ)),'k:','LineWidth',2);
end