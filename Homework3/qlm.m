function qlm = qlm(l,m,q)

%Generates Even Source Moment for l,m

if ( m == 0 )
    
    dlm = ((-1)^m)*sqrt(((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
else
    
    dlm = ((-1)^m)*sqrt(2*((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
end

plm = AssociatedLegendrePolynomial(l,m);

qlmg = @(theta,phi) dlm.*plm(cos(theta)).*cos(m.*phi).*q;

qlm = LQnQuad(qlmg,4);

return
