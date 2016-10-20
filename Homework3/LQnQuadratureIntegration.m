function I = LQnQuadratureIntegration(f,N)

[dirs,w] = LQnQuadrature(N);

sumI = 0;

for i = 1:length(w)
    
    mu = dirs(1,i);
    eta = dirs(2,i);
    xi = dirs(3,i);
    
    sumI = sumI + w(i).*f(mu,eta,xi);
    
end

I = sumI*(4*pi/8);

return