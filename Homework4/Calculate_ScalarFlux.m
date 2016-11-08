function scalar_flux = Calculate_ScalarFlux(xi,wi,angular_flux)

scalar_flux = zeros(1,length(xi));

for j = 1:length(xi)
        
    for l = 1:length(wi)
            
        scalar_flux(j) = scalar_flux(j) + (1/2)*wi(l)*angular_flux(l,j);
            
    end
        
end

return