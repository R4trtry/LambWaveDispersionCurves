% -------------------------------------------------------------------------
% dispcmaker.m
%   Generate dispersion curves by numerically solve the wave equations
%   Computation follows steps in Rose, Joseph L. Ultrasonic waves in solid
%   media. Cambridge university press, 2004.
%
% -------------------------------------------------------------------------
% Author: Chang Liu, changliu.cee@gmail.com
% Last updated: Nov 2, 2014
%
% -------------------------------------------------------------------------

clear;clc;close all
%% ------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------

% PROPERTIES YOU CAN CHANGE
% Material properties
cl = 6300;      % longitudinal wave speed
ct = 3100;      % torsional wave speed
d  = .01;        % thickness = 10 mm

% Computation properties
SAMPLING_FREQUENCY  = 2000000;       % Sampling freq = 2 MHz
FREQUENCY_STEP      = 1000;
MAX_PHASE_VELOCITY  = 20000;
PHASE_VELOCITY_STEP = .1;

% PROPERTIES, DO NOT CHANGE
% Discretize phase velocity and frequency
cp = PHASE_VELOCITY_STEP : PHASE_VELOCITY_STEP : MAX_PHASE_VELOCITY;
freq = 1 : FREQUENCY_STEP : SAMPLING_FREQUENCY/2;

Lc = length(cp);  Lf = length(freq);
fd = freq.*d;   % frequency-thickness product

% Cutoff frequency
cutoffs = [0 ct*(1:(fd(end)/ct)) cl*(.5:fd(end)/cl)];
cutoffa = [0 cl*(1:(fd(end)/cl)) ct*(.5:fd(end)/ct)];

% Initialize results
dispcs = zeros(Lf,length(cutoffs));
dispca = zeros(Lf,length(cutoffa));

%% ------------------------------------------------------------------------
% SOLVE WAVE EQUATIONS AT EACH FREQUENCIES
% -------------------------------------------------------------------------

eps2 = eps^2; t0 = tic; ts = 0;
for j = 1:Lf

    fprintf([ repmat('\b', 1, 41) '%08i / %08i [Time left: %s]'],...
              j, Lf, datestr(ts/24/3600/(j)*(Lf-j+1), 'HH:MM:SS'));  
    
    f = freq(j);
    
    % residue vector for symmetric and asymmetric modes
    sym = lamb1(f,cp,d,cl,ct);      asm = lamb2(f,cp,d,cl,ct);
    
    % in-line function of c(velocity) for further root finding
    Sym = @(c)lamb1(f,c,d,cl,ct);   Asm = @(c)lamb2(f,c,d,cl,ct);
    
    cps = []; cpa = [];
    sym(abs(sym)>0.15) = nan;        asm(abs(asm)>100) = nan;
    % check sign change
    id1 = (sym(1:end-1) .* sym(2:end) < -eps2);
    id2 = (asm(1:end-1) .* asm(2:end) < -eps2);
    
    % find the exact root in that grid, record the root as velocity
    for i = find(id1)
        czero = fzero(Sym, [cp(i),cp(i+1)]);
        cps(end+1) = czero;
    end
    
    % do the same for asym mode
    for i = find(id2)
        czero = fzero(Asm, [cp(i),cp(i+1)]);
        cpa(end+1) = czero;
    end
    
    %
    scount = length(cps);           acount = length(cpa);
    dispcs(j,1:scount) = cps;       dispca(j,1:acount) = cpa;
    
    ts = toc(t0);
end

% CHECK ZEROS / NANS / DISCONTINUITIES
dispca = correction(dispca);
dispcs = correction(dispcs);
% CONCATENATE TWO SETS OF CURVES 
dispc = [dispca dispcs];

%% ------------------------------------------------------------------------
% PLOT DISPERSION CURVES, SAVE IN MAT
% -------------------------------------------------------------------------
save('../velocity.mat','dispc','fd');
figure;hold on
plot(fd,dispc(:,1:acount),'--');    plot(fd,dispc(:,1+acount:end),'-');

