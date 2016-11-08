function [xi,scalar_flux] = OneDDiscreteOrdinates(mu,wi,h,alpha,L,sigt,sigs,qex)

xii = 0:h:L;
xi = h/2:h:L;

N = length(xi);
    
half_angular_flux = zeros(length(mu),length(xii));
angular_flux = zeros(length(mu),length(xi));

for iter = 1:1000
    
    scalar_flux_old = Calculate_ScalarFlux(xi,wi,angular_flux);
    
    for j = 1:length(mu)
        
        if ( mu(j) > 0 )
            
            %Set boundary condition on angular flux
            half_angular_flux(j,1) = 2.0;
            
            for k = 1:length(xi)
                
                qscat = 0;
                
                for l = 1:length(wi)
                    
                    qscat = qscat + (1/2)*sigs*wi(l)*angular_flux(l,k);
                    
                end
                
                q = qscat + qex;
                
                angular_flux(j,k) = (2*abs(mu(j))/(1+alpha) + sigt*h)^(-1)*(h*q + ...
                    abs(mu(j))*half_angular_flux(j,k)*(1 + (1-alpha)/(1+alpha)));
                
                half_angular_flux(j,k+1) = (2/(1+alpha))*angular_flux(j,k) - ...
                    ((1-alpha)/(1+alpha))*half_angular_flux(j,k);
                
            end
            
        elseif ( mu(j) < 0 )
            
            %Set boundary condition on angular flux
            half_angular_flux(j,N+1) = half_angular_flux(length(wi)+1-j,N+1);
            
            for k = length(xi):-1:1
                
                qscat = 0;
                
                for l = 1:length(wi)
                    
                    qscat = qscat + (1/2)*sigs*wi(l)*angular_flux(l,k);
                    
                end
                
                q = qscat + qex;
                
                angular_flux(j,k) = (2*abs(mu(j))/(1+alpha) + sigt*h)^(-1)*(h*q + ...
                    abs(mu(j))*half_angular_flux(j,k+1)*(1 + (1-alpha)/(1+alpha)));
                
                half_angular_flux(j,k) = (2/(1-alpha))*angular_flux(j,k) - ...
                    ((1+alpha)/(1-alpha))*half_angular_flux(j,k+1);
                
            end
            
        end
        
    end
    
    scalar_flux = Calculate_ScalarFlux(xi,wi,angular_flux);
    
    norm_flux = sqrt(sum((scalar_flux - scalar_flux_old).^2))
    
    if ( norm_flux < 1e-3)
        
        break
        
    end
    
end

return