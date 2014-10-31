function [asm] = lamb2(f,cp,d,cl,ct)
    w = 2*pi*f;
    h = d/2;
    
    k2 = (w./cp).^2;
    p = sqrt((w/cl)^2 - k2);
    q = sqrt((w/ct)^2 - k2);
    
    tph = tan(p.*h);
    tqh = tan(q.*h);
    qk = (q.^2-k2).^2;
    
    asm = (tqh.*q + qk.*tph./4./k2./p);
end