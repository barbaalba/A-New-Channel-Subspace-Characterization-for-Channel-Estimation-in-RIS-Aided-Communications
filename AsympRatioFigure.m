% To generate the asymp ratio figure 
clear; clc; close all;
azref = pi/2; % Initial point
elref = 0; % Initial point
M_V = transpose(4*2.^(0:5)); M_H = M_V; M = M_H.*M_V; % Antennas config
d_V = [1/2; 1/4]; d_H = d_V; % Antennas config
% collect the DOF according to the proposed algorithm
bsnum = zeros(length(M_V),length(d_V));
for j = 1:length(d_V)
    for i = 1:length(M_H)
        [~,~,bsnum(i,j)] = UPA_BasisElup(M_V(i),M_H(i),d_V(j),d_H(j),azref,elref);
    end
end
Ratio = bsnum./M;
AsympRatio = pi*d_H.*d_V;
AsympRatio = repelem(AsympRatio,size(Ratio,1));
AsympRatio = reshape(AsympRatio,length(M_V),length(d_V));
figure('DefaultAxesFontSize',20,'defaultLineLineWidth',2,'defaultAxesTickLabelInterpreter','latex');
for i = 1:length(d_V)
    plot(Ratio(:,i),'-o','MarkerSize',10);
    hold on;
    plot(AsympRatio(:,i),'--','marker','diamond','MarkerSize',10);
end
grid on;
ylabel('$\frac{\eta}{M}$','Interpreter','latex','FontSize',20);
ylim([0,1]);
yticks(0:0.1:1);
xlim([1,length(M_H)]);
xticks(1:6);
xticklabels({'$4 \times 4$','$8 \times 8$','$16 \times 16$','$32 \times 32$','$64 \times 64$','$128 \times 128$'});
legend('Proposed, $\lambda/2$','Asymp, $\lambda/2$','Proposed, $\lambda/4$','Asymp, $\lambda/4$','interpreter','latex');
%% Figure Configuration
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);