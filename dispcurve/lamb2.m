% -------------------------------------------------------------------------
% lamb2.m
%   lamb wave equation (residual) for asymetric mode. 
%   Computation follows steps in Rose, Joseph L. Ultrasonic waves in solid
%   media. Cambridge university press, 2004.
%
% -------------------------------------------------------------------------
% Author: Chang Liu, changliu.cee@gmail.com
% Last updated: Nov 2, 2014
%
% ------------------------------------------------------------------------- 
% INPUT:
%   -f:  frequency, discretized [Hz]
%   -cp: phase velocity, discretized [m/s]
%   -d:  thickness [m]
%   -cl: longitudinal wave velocity [m/s]
%   -ct: tortional wave velocity [m/s]
%
% OUTPUT:
%   -asm: function handle for the residual of wave equation
% -------------------------------------------------------------------------

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