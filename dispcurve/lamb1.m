function [sym] = lamb1(f,cp,d,cl,ct)
    w = 2*pi*f;
    h = d/2;
    
    k2 = (w./cp).^2;
    p = (sqrt((w/cl)^2 - k2));
    q = (sqrt((w/ct)^2 - k2));
    
    tph = tan(p.*h);
    tqh = tan(q.*h);
    qk = (q.^2-k2).^2;
    
    sym = (tqh./q + 4.*k2.*p.*tph./qk);
end