function W = DFTBookBuild(M_H,M_V)
% This function produces DFT code book for a UPA antenna using Kronecker product
% with any antenna spacing
% M_H: number of antenna elements on horizental axis
% M_V: number of antenna elements on vertical axis
% W: DFT precoder

y = zeros(M_H,M_H); % DFT for Horizontal ULA
z = zeros(M_V,M_V); % DFT for Vertical ULA

% Horizental DFT
x = exp(-1i*2*pi/M_H);
for i = 1:M_H
    for j = 1:M_H
       y(i,j) = x^((i-1)*(j-1));
    end
end

% Vertical DFT
x = exp(-1i*2*pi/M_V);
for i = 1:M_V
    for j = 1:M_V
        z(i,j) = x^((i-1)*(j-1));
    end
end

% DFT Code book for UPA
W = kron(z,y);

end