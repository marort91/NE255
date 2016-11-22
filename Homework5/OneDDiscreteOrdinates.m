function [xi,scalar_flux] = OneDDiscreteOrdinates(mu,wi,h,NEgrps,alpha,L,sigt,sigs,qex,tol)

xii = 0:h:L;
xi = h/2:h:L;

N = length(xi);
    
half_angular_flux = zeros(length(mu),length(xii),NEgrps);
angular_flux = zeros(length(mu),length(xi),NEgrps);

half_angular_flux(1:length(mu)/2,1,1) = 0.5;

for iter = 1:1000
    
    %scalar_flux_old = Calculate_ScalarFlux(xi,wi,NEgrps,angular_flux);
    scalar_flux_old_total = Calculate_ScalarFlux(xi,wi,NEgrps,angular_flux);
    
    for e = 1:NEgrps
        
    for InIter = 1:1000;
        
        scalar_flux_old = Calculate_ScalarFlux(xi,wi,NEgrps,angular_flux);
    
    for j = 1:length(mu)
        
        if ( mu(j) > 0 )
            
            %Set boundary condition on angular flux
            %half_angular_flux(j,1,1) = 0.5;
            
            for k = 1:length(xi)
                
                qscat = 0;
                 
                for srcE = 1:NEgrps
                    
                for l = 1:length(wi)
                    
                    qscat = qscat + (1/2)*sigs(e,srcE)*wi(l)*angular_flux(l,k,srcE);
                    
                end
                
                end
                
                q = (1/1)*qscat + qex(e);
                
                angular_flux(j,k,e) = (2*abs(mu(j))/(1+alpha) + sigt(e)*h)^(-1)*(h*q + ...
                    abs(mu(j))*half_angular_flux(j,k,e)*(1 + (1-alpha)/(1+alpha)));
                
                half_angular_flux(j,k+1,e) = (2/(1+alpha))*angular_flux(j,k,e) - ...
                    ((1-alpha)/(1+alpha))*half_angular_flux(j,k,e);
                
            end
            
        elseif ( mu(j) < 0 )
            
            %Set boundary condition on angular flux
            half_angular_flux(j,N+1,e) = half_angular_flux(length(wi)+1-j,N+1,e);
            
            for k = length(xi):-1:1
                
                qscat = 0;
                 
                for srcE = 1:NEgrps
                    
                for l = 1:length(wi)
                    
                    qscat = qscat + (1/2)*sigs(e,srcE)*wi(l)*angular_flux(l,k,srcE);
                    
                end
                
                end
                
                q = (1/1)*qscat + qex(e);
                
                angular_flux(j,k,e) = (2*abs(mu(j))/(1+alpha) + sigt(e)*h)^(-1)*(h*q + ...
                    abs(mu(j))*half_angular_flux(j,k+1,e)*(1 + (1-alpha)/(1+alpha)));
                
                half_angular_flux(j,k,e) = (2/(1-alpha))*angular_flux(j,k,e) - ...
                    ((1+alpha)/(1-alpha))*half_angular_flux(j,k+1,e);
                
            end
            
        end
        
    end
    
    scalar_flux = Calculate_ScalarFlux(xi,wi,NEgrps,angular_flux);
    
    norm_flux = sqrt(sum((scalar_flux(e,:) - scalar_flux_old(e,:)).^2));
    
    fprintf('Energy Group %i, Iteration: %i, Norm: %f \n',e,InIter,norm_flux);
    
    if ( norm_flux < tol)
        
        break
        
    end
    
    end
    
    end
    
    norm_flux_total = sqrt(sum(scalar_flux - scalar_flux_old_total).^2);
    
     if ( norm_flux_total < tol)
        
        break
        
    end
    
end

return