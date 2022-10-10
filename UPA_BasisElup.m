function [theta,phi,bsnum] = UPA_BasisElup(M_V,M_H,d_V,d_H,azref)
% This function evaluate the orthogonal basis for UPA structure
% Input
%   M_V : number of antennas in vertical axis
%   M_H : number of antennas in horizental axis
%   d_V : normalized antenna spacing in vertical axis
%   d_H : normalized antenna spacing in horizental axis
%   azref : the reference azimuth angle
% Output
%   theta,phi : angle pairs
%   bsnum : number of orthogonal beams

L_V = M_V * d_V;
L_H = M_H * d_H;

%% first find the elevation orthogonal to 90 degree
theta = [0]; % 0 as the elevation reference
for k = 1:M_V - 1 
    omega = k / L_V;
    if mod(omega,1/d_V) <= 1
        theta = [theta asin(omega)];
    elseif mod(omega,-1/d_V) >= -1
        omega = mod(omega,-1/d_V);
        theta = [theta asin(omega)];
    end
end

%% Initialize the azimuth angles with the reference angle
phi = struct;
for j = 1:length(theta)
    eval(['phi.theta' num2str(j) '= [azref];']);
end
fn = fieldnames(phi);

%% Build the azimuth for each elevation angle
for j = 1:length(theta)

    for k = 1:M_H-1

        gama = k/L_H;

        if mod(gama,1/d_H) <= 2
            gama = gama/cos(theta(j)) + sin(azref);
            val = asin(gama);
            if isreal(val)
                phi.(fn{j}) = [phi.(fn{j}) val];
            end
        end
        gama = k/L_H;
        if mod(gama,-1/d_H) >= -2
            gama = mod(gama,-1/d_H)/cos(theta(j)) + sin(azref);
            val = asin(gama);
            if isreal(val)
                phi.(fn{j}) = [phi.(fn{j}) val];
            end
        end
        
    end

end
fn = fieldnames(phi);
bsnum = 0; % Basis number
for j = 1:length(fn)
    bsnum = bsnum + length(phi.(fn{j}));
end
disp(['Number of basis before process: ' num2str(bsnum)]);

end