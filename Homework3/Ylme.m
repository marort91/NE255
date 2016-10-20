function Ylme = Ylme(l,m)

%Generates Even Spherical Harmonic for l,m

if ( m == 0 )
    
    dlm = ((-1)^m)*sqrt(((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
else
    
    dlm = ((-1)^m)*sqrt(2*((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
end

plm = AssociatedLegendrePolynomial(l,m);

Ylme = @(theta,phi) dlm.*plm(cos(theta)).*cos(m.*phi);

return