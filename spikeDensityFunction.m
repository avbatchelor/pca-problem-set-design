function [sdf, rasterSum] = spikeDensityFunction(raster,binSize,sampleRate,sigma)
% raster: MxN where M=repeats, N=time

rasterSum = sum(raster,1);
binSamples = ceil((binSize / 1000) * sampleRate);
gaussWidth = binSamples;

if exist('sigma','var')
    sigma = ceil((sigma / 1000) * sampleRate);
else
    sigma = binSamples;
end

gSamp = -round(3*gaussWidth/2):round(3*gaussWidth/2);
% kern = normpdf(gSamp,0,sigma);
kern = 1/(sigma * sqrt(2*pi)).*exp(-gSamp.^2/(2*sigma.^2));
sdf = conv(rasterSum,kern);
sdf = sdf(1:length(rasterSum));