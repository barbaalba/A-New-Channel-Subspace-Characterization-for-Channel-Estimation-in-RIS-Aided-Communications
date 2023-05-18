%% This code is written to implement the end-end beam construction with RIS 
% using the basis vectors generated from physical structure 
% Assumption: BS and RIS are located in the same elevation angle right in
% front of each other
clear; clc; close all;
%% Initialization 
Ptx = 20; %[dBm]
x = 1*sqrt(db2pow(Ptx-30)); % Tx signal
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_V = 128; % Vertical elements of RIS
M_H = 128; % Horizontal elements of RIS
M = M_V * M_H; % total number of antennas in URA
RISspacing = 1/4; % antennas spacing with respect to lambda RIS
numUE = 2;
N = 64; % Antenna numbers of BS
BSspacing = 1/2;
noisepowdB = -96-30;
noisepow = db2pow(noisepowdB); 
% channel parameters (LOS mmWave) alpha + 10*beta*log10(d)
alpha = 61.4;
beta = 2;
plEval = @(d) alpha + 10*beta*log10(d);

%% Generate Channel between RIS and UE
RIS_AOD = unifrnd(-pi/2,pi/2,1,numUE); % azimuth of departure from RIS to UE
RIS_EOD = unifrnd(-pi/2,0,1,numUE); % Elevation of departure from RIS to UE
d = unifrnd(30,50,1,numUE); % distance between user and RIS
UEPLdB = plEval(d); % Evaluate the corresponding attenuation
phshift = exp(-1i*2*pi/lambda*d); % the corresponding phase shift due to propagation
Chancoeff = sqrt(db2pow(-UEPLdB)) .* phshift; % the channel complex coefficient

h = Chancoeff.*UPA_Evaluate(lambda,M_V,M_H,RIS_AOD,RIS_EOD,RISspacing,RISspacing);; % Uplink Channel (RIS-UE)

%Add NLOS components
Ricefactor = 13 - 0.03*d; % in dB
parfor n = 1:numUE
    NLOSp = Ricianpow(db2pow(Ricefactor(n)),db2pow(-UEPLdB(n)));
    Chancoeff = sqrt(NLOSp/2) .* (randn(1,length(NLOSp)) + 1i*randn(1,length(NLOSp)));
    RIS_AOD = unifrnd(-pi/2,pi/2,1,length(NLOSp));  
    RIS_EOD = unifrnd(-pi/2,0,1,length(NLOSp)); 
    A = UPA_Evaluate(lambda,M_V,M_H,RIS_AOD,RIS_EOD,RISspacing,RISspacing);
    h(:,n) = h(:,n) + sum(Chancoeff.*A,2);
end
%% Generate channel between RIS and BS

BS_AOD = 0;%unifrnd(-pi,pi/2,1); % Azimuth of departure from BS to RIS

RIS_AOA = pi/6; % Azimuth of arrival at RIS from BS
RIS_EOA = 0; % Elevation of arrival to RIS from BS (They are in the same level)
Wbs = 1/sqrt(N) * ULA_Evaluate(lambda,N,BS_AOD,BSspacing); % norm-1
A = ULA_Evaluate(lambda,N,BS_AOD,BSspacing) * ...
    UPA_Evaluate(lambda,M_V,M_H,RIS_AOA,RIS_EOA,RISspacing,RISspacing)';
d = 10;%unifrnd(1,10,1);
PLdB = plEval(d); % channel attenuation
phshift = exp(-1i*2*pi/lambda*d); % channel phase shift due to propagation
ChanComp = sqrt(db2pow(-PLdB)) * phshift; % channel complex coefficient

H =  ChanComp * A; % Uplink Channel (BS-RIS)

% Build the cascaded channel
V =  zeros(M*N,numUE);
parfor j = 1:numUE
    Vred = H*diag(h(:,j)); % redundancy variable
    V(:,j) = reshape(Vred,[],1);
end
%% Establish the basis angles
[theta,phi,bsnum] = UPA_BasisElupnew(M_V,M_H,RISspacing,RISspacing,pi/2,0);
fn = fieldnames(phi);
W = UPA_Codebook(lambda,theta,phi,M_V,M_H,RISspacing,RISspacing); % Set of basis vector
W = W.*UPA_Evaluate(lambda,M_V,M_H,RIS_AOA,RIS_EOA,RISspacing,RISspacing);
%% Chanenl estimation through our algorithm 

% Sensing UEs
y = zeros(size(W,2),numUE);
parfor i = 1:size(W,2)
    noise = sqrt(noisepow/2) * (randn(N,numUE) + 1i * randn(N,numUE));
    y(i,:) = Wbs' * (H * diag(W(:,i)) * h * x + noise);
end
% estimate
Wopt = W * conj(y) / M / sqrt(db2pow(Ptx-30)); % Eq. 27
Vest =  zeros(M*N,numUE);
parfor i = 1:numUE
    Vest(:,i) = reshape(Wbs * Wopt(:,i)',[],1); % Eq. 28
end
Wopt = angle(Wopt);
Wopt = exp(1i*Wopt); % Eq. 29
NMSE_V = mean(vecnorm(V-Vest).^2./vecnorm(V).^2); % Eq. 30
disp(['Porposed estimation NMSE: ' num2str(pow2db(NMSE_V))]);
%% LS estimator using DFT config for RIS
Wd = BuildDFT(M); % M*M matrix used for training (row per time)
Psi = kron(Wd,eye(N)); 
X = diag(repelem(x,N*M));
Comb = X*Psi;
% Estimate the cascaded channel
rx = Comb * V + sqrt(noisepow/2) * (randn(N*M,numUE) + 1i * randn(N*M,numUE));
Vest = (Comb'*Comb)\(Comb') * rx; % Estimated cascaded channel 
LSNMSE = mean(vecnorm(V-Vest).^2./vecnorm(V).^2);
disp(['End to End LS NMSE: ' num2str(pow2db(LSNMSE))]);
parfor j = 1:numUE
    Vred = Vest(:,j);
    Vred = reshape(Vred,N,M);
    WLS(:,j) = angle(Vred'*Wbs);
    WLS(:,j) = exp(1i*WLS(:,j)); % Final RIS config using LS estimator
 end
%% Measure the performance
% Measure rx power using the configurations and best config (Known Channel)
rx = zeros(1,numUE); % Collect Rx for our method
rxbest = zeros(1,numUE); % Collect Rx for optimum config
rxDFT = zeros(numUE,M); % Collect Rx for LS 
coeff = zeros(M,numUE); % Optimum RIs configuration
parfor j = 1:numUE
    V = H*diag(h(:,j));
    thet = angle(Wbs'*V);
    coeff(:,j) = exp(-1i*thet); % optimum RIS config matching cascaded channel
    rx(j) = abs(Wbs' * H * diag(Wopt(:,j)) * h(:,j) * x).^2; % our method
    rxbest(j) = abs(Wbs' * H * diag(coeff(:,j)) * h(:,j) * x).^2; % Optimum
    rxLS(j) = abs(Wbs' * H * diag(WLS(:,j)) * h(:,j) * x).^2; % LS 
end

% Evaluate the SNR
rxdB = pow2db(rx);
SNRdB = mean(rxdB) - noisepowdB;
disp(['The measured SNR is:' num2str(SNRdB)]);

rxdB = pow2db(rxbest);
SNRdB = mean(rxdB) - noisepowdB;
disp(['The max SNR is:' num2str(SNRdB)]);

rxLSdB = pow2db(rxLS);
SNRdB = mean(rxLSdB) - noisepowdB;
disp(['The Least Square measured SNR is:' num2str(SNRdB)]);
