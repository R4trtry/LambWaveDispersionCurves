clear;clc;close all

%% dispTry2
timer = tic;

%%
% ------------------------------------------------------------------------
% PARAMETERS
% ------------------------------------------------------------------------

% Material properties
cl = 6300;      % longitudinal wave speed
ct = 3100;      % torsional wave speed
d = .01;        % thickness = 10 mm

% Computation properties
SAMPLING_FREQUENCY  = 4000000;       % Sampling freq = 2 MHz
FREQUENCY_STEP      = 1000;
MAX_PHASE_VELOCITY  = 20000;
PHASE_VELOCITY_STEP = .1;

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

%%
% ------------------------------------------------------------------------
% SOLVE WAVE EQUATIONS AT EACH FREQUENCIES
% ------------------------------------------------------------------------

h = waitbar(0,'Progress','Name','Imaging',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

eps2 = eps^2;t2 = tic;
for j = 1:Lf
    if getappdata(h,'canceling'), break; end
    waitbar(j/Lf,h,sprintf('Progress %.2f%%, Remain time = %.1f sec',...
        j/Lf*100,toc(t2)*(Lf/j-1)))
    
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
end
delete(h)

%%
% ------------------------------------------------------------------------
% CHECK ZEROS / NANS / DISCONTINUITIES
% ------------------------------------------------------------------------
dispca(dispca<dispca(2,1))=NaN;dispca(1,1) = 0;
dispcs(dispcs<dispca(2,1))=NaN;
correction; % fix the curve by sorting and matching the order
dispc = [dispca dispcs];

%% 
% ------------------------------------------------------------------------
% PLOT DISPERSION CURVES, SAVE IN MAT
% ------------------------------------------------------------------------

toc(timer)
save('../velocity.mat','dispc','fd');
figure;hold on
plot(fd,dispc(:,1:acount),'--');    plot(fd,dispc(:,1+acount:end),'-');