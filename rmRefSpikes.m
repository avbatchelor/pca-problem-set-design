function [raster,didShift] = rmRefSpikes(raster,refPerBins,shiftFlag)

for i = refPerBins+1:numel(raster)
    % If there are spikes in the refractory period bins set them to zero
    if sum(raster((i-refPerBins):i)) > 1
        raster((i-refPerBins+1):i) = 0;
    end
end

% Shift by a single bin to avoid periodicity in modeled psd estimates
if shiftFlag
    didShift = randi([-1 1]);
    raster = circshift(raster,[0,didShift]);
else
    didShift = 0;
end