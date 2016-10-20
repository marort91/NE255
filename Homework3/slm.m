function slm = slm(l,m,q)

%Generates Odd Source Moment for l,m

if ( m == 0 )
    
    dlm = ((-1)^m)*sqrt(((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
else
    
    dlm = ((-1)^m)*sqrt(2*((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
end

plm = AssociatedLegendrePolynomial(l,m);

slmg = @(theta,phi) dlm.*plm(cos(theta)).*sin(m.*phi).*q;

slm = LQnQuad(slmg,4);

return