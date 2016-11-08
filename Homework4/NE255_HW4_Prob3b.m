%Problem 3
%NE 255
%1D Discrete Ordinates Code

clc, clear, clf

L = 2.0;

alpha = [-0.5 0 0.5];
mu = [0.7 0.2 -0.2 -0.7];
wi = [0.5 0.5 0.5 0.5];
sigt = 1.0;
sigs = 0.9;
qex = 1.0;

h = 0.08;

for i = 1:length(alpha)

    [xi,scalar_flux] = OneDDiscreteOrdinates(mu,wi,h,alpha(i),L,sigt,sigs,qex);
    plot(xi,scalar_flux)
    grid on
    xlabel('Distance (cm)');
    ylabel('Scalar Flux');
    titl = sprintf('Scalar Flux for Sigs = %.2f, Q = %.2f, Alpha = %.2f',sigs,qex,alpha(i));
    title(titl);
    fid = sprintf('scattering_alpha%i',i);
    export_fig(fid,'-pdf','-nocrop')
    
end