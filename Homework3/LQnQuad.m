function I = LQnQuad(f,N)

[dirs,w] = LQnQuadrature(N);

mu = dirs(1,:);
eta = dirs(2,:);
xi = dirs(3,:);

theta = acos(mu);

for i = 1:length(xi)
    
    phi(i) = atan2(xi(i),eta(i));
    
end

sumI = 0;

for i = 1:length(w)
    
    sumI = sumI + w(i)*f(theta(i),phi(i));
    
end

I = sumI*(4*pi/8);

return