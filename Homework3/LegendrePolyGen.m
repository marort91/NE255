function pn = LegendrePolyGen(N)

p0 = @(x) 1;
p1 = @(x) x;

if ( N == 0 )
    
    pn = p0;
    
elseif ( N == 1 )
    
    pn = p1;
    
else

    for i = 1:N-1
    
        p2 = @(x) (2*i+1).*x.*p1(x) - i.*p0(x);
        p2 = @(x) p2(x)./(i+1);
    
        p0 = p1;
        p1 = p2;
    
    end

    pn = p2;
    
end

return