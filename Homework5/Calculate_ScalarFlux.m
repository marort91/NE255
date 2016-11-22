function scalar_flux = Calculate_ScalarFlux(xi,wi,Negrp,angular_flux)

scalar_flux = zeros(Negrp,length(xi));

for e = 1:Negrp

for j = 1:length(xi)
        
    for l = 1:length(wi)
            
        scalar_flux(e,j) = scalar_flux(e,j) + (1/2)*wi(l)*angular_flux(l,j,e);
            
    end
        
end

end

return