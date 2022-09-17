function [theta,phi] = UPA_BasisEl(lambda,M_V,M_H,d_V,d_H)
% This function evaluate the orthogonal basis for UPA structure
L_V = M_V * d_V;
L_H = M_H * d_H;

%% first find the elevation orthogonal to 90 degree
theta = [0];
for k = 1:M_V - 1 

    omega = k / L_V;

    if mod(omega,1/d_V) <= 1

        theta = [theta asin(omega)];

    elseif mod(omega,-1/d_V) >= -1

        theta = [theta asin(mod(omega,-1/d_V))];

    end

end

%% Initialize the azimuth angles
phi = struct;
for j = 1:length(theta)

    eval(['phi.theta' num2str(j) '= [0];']);

end
fn = fieldnames(phi);

%% Build the azimuth for each elevation angle
for j = 1:length(theta)

    for k = 1:M_H-1

        gama = k/L_H;

        if mod(gama,1/d_H) <= 1
            gama = gama/cos(theta(j));
            val = asin(gama);
            if isreal(val)
                phi.(fn{j}) = [phi.(fn{j}) val];
            end
        elseif mod(gama,-1/d_H) >= -1
            gama = mod(gama,-1/d_H)/cos(theta(j));
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