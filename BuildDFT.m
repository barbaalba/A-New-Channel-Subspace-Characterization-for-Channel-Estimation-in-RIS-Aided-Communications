function o = BuildDFT(M)
% To build a DFT precoder matrix according to 
% T. L. Jensen and E. De Carvalho, "An Optimal Channel Estimation Scheme 
% for Intelligent Reflecting Surfaces Based on a Minimum Variance Unbiased 
% Estimator," ICASSP 2020 - 2020 IEEE International Conference on 
% Acoustics, Speech and Signal Processing (ICASSP), 2020
o = zeros(M,M);
for t = 1:M
    for m = 1:M
        o(t,m) = exp(-1i*2*pi/M*(t-1)*(m-1));
    end   
end
end