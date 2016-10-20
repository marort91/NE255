function [dirs,wi] = LQnQuadrature(SN)

if ( SN == 4 )
    
    mu(1) = 0.3500212;
    mu(2) = 0.8688903;
    mu(3) = -mu(2);
    mu(4) = -mu(1);
    
    eta = mu;
    xi = mu;
    
    ordinates = combvec(mu,eta,xi);
    mag = round(sum(ordinates.^2),6);

    locs = mag == 1;

    dirs = ordinates(:,locs);
    wi = 0.3333333.*ones(1,length(dirs(1,:)));

elseif ( SN == 6 )
    
    mu(1) = 0.2666355;
    mu(2) = 0.6815076;
    mu(3) = 0.9261808;
    mu(4) = -mu(3);
    mu(5) = -mu(2);
    mu(6) = -mu(1);
    
    eta = mu;
    xi = mu;
    
    w(1) = 0.1761263;
    w(2) = 0.1572071;
    
    ordinates = combvec(mu,eta,xi);
    mag = round(sum(ordinates.^2),6);

    locs = mag == 1;

    dirs = ordinates(:,locs);
    
    for i = 1:length(dirs)
        
        if ( abs(dirs(3,i)) == max(xi) )
            
            dirs(4,i) = w(1);
            
        elseif ( abs(dirs(3,i)) == median(abs(xi)) )
            
            dirs(4,i) = w(2);
            
        elseif ( abs(dirs(3,i)) == min(abs(xi)) && abs(dirs(2,i)) == median(abs(eta)) &&...
                abs(dirs(1,i)) == median(abs(mu)) )
            
            dirs(4,i) = w(2);
            
        else
            
            dirs(4,i) = w(1);
            
        end
        
    end
    
    wi = dirs(4,:);
    dirs = dirs(1:3,:);
    
else
    
    print('Higher Order LQn Sets to be Implemented')
    
end


return

