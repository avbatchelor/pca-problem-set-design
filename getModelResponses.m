function [rasters,pools] = getModelResponses(filts,weights,stim,baseFR,modFR,binSize,nReps,sampleRate)
% return rasters cell array and pools matrix
[M,N] = size(stim);
pools = zeros([M N]);

% Set firing rates (this means modFR needs to be > baseFR by at least 20%)
baseFR  = baseFR + 0.1*baseFR*rand(1);
modFR   = modFR + 0.1*modFR*rand(1) - baseFR;

for iS = 1:M
    for iF = 1:numel(filts)
        
        % Convolve all filters (unused ones have zero weight)
        output = conv(stim(iS,:),filts{iF});
        output = output(1:N)./max(output(1:N));
        
        % Rectify output
        output(output < 0) = 0;
        
        pools(iS,:) = pools(iS,:) + (output .* weights{iF});
    end
    
    % Rectify pool
    pools(iS,:) = pools(iS,:) .* (pools(iS,:) > 0);
    
    % Downsample pool to binsize... probably a crappy way to do this
    binSamps = (binSize/1000) * sampleRate;
    nBins = length(pools(iS,:))/binSamps;
    dsPool = downsample(pools(iS,:),binSamps);
    
    % Offset pool to approximate baseline firing rate
    dsPool = dsPool + (baseFR/sampleRate)*binSize;
    
    % Poisson (refractory corrected) spiking output
    rasters{iS} = zeros(nReps,nBins);
    for iR = 1:nReps
        rasters{iS}(iR,:) = rand(1, nBins) < dsPool.*(modFR);
        rasters{iS}(iR,:) = rmRefSpikes(rasters{iS}(iR,:),1,1);
    end
end
