function p = Ricianpow(Rfact,LOSpow)
lambda = 1.8; % For Cluster 
rdelay = 2.8; % Constant variable 
zeta = 4; % Constant variable

NLOSpow = LOSpow / Rfact;
K = max(1,poissrnd(lambda)); % Number of clusters
gamap = exp((1- rdelay) * rand(1,K));
gamap = gamap .* 10.^(-0.1 * zeta * randn(1,K));
gama = sort(gamap / sum(gamap),'descend');
p = NLOSpow * gama;
end