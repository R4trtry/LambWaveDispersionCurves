% -------------------------------------------------------------------------
% correction.m
%   check and fix ZEROS / NANS / DISCONTINUITIES
%   update dispc in place
%
% -------------------------------------------------------------------------
% Author: Chang Liu, changliu.cee@gmail.com
% Last updated: Nov 2, 2014
%
% -------------------------------------------------------------------------

function dispc = correction(dispc)
    dispc(dispc<eps)=NaN;
    [numFreq, numCurves] = size(dispc);
    
    % loop through frequencies
    for j = numFreq-1 : -1 : 1
        
        % loop through curves
        for i = numCurves: -1 : 1
            
            % 
            [~,I]  = min(abs(dispc(j  ,i) - dispc(j+1,:)));
            [~,I1] = min(abs(dispc(j+1,I) - dispc(j  ,:)));

            if(~isnan(I) && I~=i && I1==i), 
                dispc(j,I) = dispc(j,i); dispc(j,i) = nan; 
            end

        end

    end

    
end