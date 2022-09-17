function o = BuildDFT(M)
o = zeros(M,M);
for t = 1:M
    for m = 1:M
        o(t,m) = exp(-1i*2*pi/M*(t-1)*(m-1));
    end   
end
end