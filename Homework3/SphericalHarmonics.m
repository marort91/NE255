function [ylm,f] = SphericalHarmonics(l,m,meshN)

clf

if ( m < 0 )
    
    T0 = ((-1)^(-m))*(factorial(l-(-m))/factorial(l+(-m)));
    
else
    
    T0 = 1;
    
end

[theta,phi] = meshgrid(linspace(0,pi,meshN),linspace(0,2*pi,meshN));

plm = AssociatedLegendrePolynomial(l,m);

%Generate Even and Odd Spherical Harmonic Functions

if ( m == 0 )
    
    dlm = ((-1)^m)*sqrt(((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
else
    
    dlm = ((-1)^m)*sqrt(2*((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));
    
end

ylme = @(theta,phi) dlm.*plm(cos(theta)).*cos(m.*phi);
ylmo = @(theta,phi) dlm.*plm(cos(theta)).*sin(m.*phi);

T1 = (-1)^m;
T2 = sqrt( ((2*l+1)/(4*pi))*(factorial(l-m)/factorial(l+m)));

f = @(theta,phi) T0.*T1*T2*plm(cos(theta)).*exp(1i.*m.*phi);

ylm = abs(f(theta,phi));
x = ylm.*sin(theta).*cos(phi);
y = ylm.*sin(theta).*sin(phi);
z = ylm.*cos(theta);

surf(x,y,z);
xlabel('x')
ylabel('y')
str = sprintf('$$|Y_{%i}^{%i}|$$',l,m);
title(str,'Interpreter','Latex','FontSize',16);
axis equal

return