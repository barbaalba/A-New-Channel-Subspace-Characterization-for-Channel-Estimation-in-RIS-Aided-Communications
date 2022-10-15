function beamforming3Dplot(lambda,M_H,M_V,d_H,d_V,Azimuth,Elevation)
% 
%Prepare to plot colors on a sphere
N = 500;
[X,Y,Z] = sphere(N);
%Prepare to compute channel gains on the sphere
gainMap = zeros(size(X));
%Go through all azimuth and elevation angles
for n = 1:size(X,1)
    for m = 1:size(X,2)
    
        %Compute received power according to (7.28) in "Massive MIMO networks"
        [phi2,theta2] = cart2sph(X(n,m),Y(n,m),Z(n,m));
        gainMap(n,m) = abs(1/sqrt(M_H*M_V)*...
            UPA_Evaluate(lambda,M_V,M_H,phi2,theta2,d_V,d_H)'*...
            UPA_Evaluate(lambda,M_V, M_H,Azimuth,Elevation,d_V,d_H)).^2;
   
    end
end

% Plot beamforming gain
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
figure;
surf(X,Y,Z,pow2db(gainMap),'EdgeColor','none');
xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
caxis([-20 pow2db(M_H*M_V)]); % limit the color map
colormap(flipud(hot));
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');
axis equal;
set(gca,'color',[0.9 0.9 0.9]);
set(gca,'fontsize',18);
view(122,30);
hold on;
% Plot the dot lines around the sphere
varphiAngles = linspace(-pi,pi,100);
x_circ = cos(varphiAngles);
y_circ = sin(varphiAngles);
plot3(x_circ,y_circ,zeros(size(x_circ)),'k:','LineWidth',2);
end
