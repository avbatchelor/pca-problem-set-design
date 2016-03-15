%%
sampleRate = 5000;
fltLen2 = ((101)/1000)*sampleRate;
fltLen = ((50)/1000)*sampleRate;
hwDg = 2;                          % Half-width of filter
tFilt = linspace(-hwDg,hwDg,fltLen2); % filter timebase
stdGaus = 0.75;                       % std of Gaussian, in deg.
gauss = exp(-tFilt.^2/stdGaus.^2);   % Gaussian envelope
mHatF = diff(diff(gauss));
mHatF = downsample(mHatF,2);
mHatF = mHatF(1:fltLen);
mHatF = mHatF./-(min(mHatF));
alphaF= alphaF./abs(max(alphaF));
expF = (1/5).^linspace(0,1,fltLen);
combF = lineF .* mHatF;

close all
plot(lineF)
hold on
plot(mHatF)
plot(combF)
legend('alpha','mhat','comb')
%%
close all
nx=80;              %Number of spatial samples in the filter
max_x =2.0;         %Half-width of filter (deg)
dx = (max_x*2)/nx;  %Spatial sampling interval of filter (deg)
 
% A row vector holding spatial sampling intervals
x_filt=linspace(-max_x,max_x,nx);
 
% Spatial filter parameters
sx=0.5;   %standard deviation of Gaussian, in deg.
sf=1.1;  %spatial frequency of carrier, in cpd
 
% Spatial filter response
gauss=exp(-x_filt.^2/sx.^2);          %Gaussian envelope
sw = gauss.*x_filt+.5;
even_x=cos(2*pi*sf*x_filt).*gauss;   %Even Gabor
odd_x=sin(2*pi*sf*x_filt).*gauss;    %Odd Gabor

plot(sw)
%plot(even_x.*x_filt)
%%
%% 2 - negative of differentiating filter with delay
output = conv(currStim,[zeros(1,20), -1 1],'same');
output = output/max(output);
% rectify output
output(output < 0) = 0;
sub(2).output = output;
sub(2).weight = (-1)^randi([0 1]) * rand(1) * .75;

% 3 - Moving average-type filter
windowSize = 26;
hannWind = hann(windowSize);
output = conv(currStim,[zeros(1,windowSize),hannWind'],'same');
output = output/max(output);
% rectify output
output(output < 0) = 0;
sub(3).output = output;
sub(3).weight = (-1)^randi([0 1]) * rand(1) * .125;

%%
sampleRate = 5000;
% stim[n=ms]
stimCmd = [0.*ones(300,1); ones(500,1); 0.*ones(300,1)];
% lowpass bessel filter the stimulus
[z,p,k] = besself(5,5000);                   
% Analog to digital mapping
[zd,pd,kd]=bilinear(z,p,k,sampleRate);      
[sos] = zp2sos(zd,pd,kd);           
% Convert to SOS form
stimFilt = sosfilt(sos,stimCmd);

% lpf = designfilt('lowpassfir', 'PassbandFrequency', 500,...
%                                 'StopbandFrequency', 800,...
%                                 'PassbandRipple', 1,...
%                                 'StopbandAttenuation', 60,...
%                                 'SampleRate', 5000,...
%                                 'DesignMethod', 'kaiserwin');
% stimFilt = filtfilt(lpf,stimCmd);

% simple haar wavelet
%cwtWeights = cwt(stimFilt,[1,2,4,5,6,7,8,9,10],'sym1');
cwtWeights = cwt(stimFilt,[59],'gaus1');
%cwtWeights = abs(cwt(stimFilt,[1,2,4,5,6,7,8,9,10],'db1'));
% sum and normalize weights
cwtWeights = (sum(cwtWeights,1));
cwtWeights = cwtWeights/max(cwtWeights);
cwtWeights = 100*cwtWeights + 1;
fr = 15; % Hz
dt = 1/1005; % s
nBins = length(cwtWeights);
for i = 1:100
    myPoissonSpikeTrain(i,:) = rand(1, nBins) < cwtWeights.*fr*dt;
end

figure
ax1=subplot(211);
plotSpikeRaster(myPoissonSpikeTrain,'PlotType','VertLine');
% plot(stimCmd)

ax2=subplot(212);
plot(stimFilt);
hold on;
plot(sum(cwtWeights,1))
linkaxes([ax1 ax2],'x');
% plotTrain = myPoissonSpikeTrain+0;
% plotTrain(plotTrain==0) = nan;
% plot(plotTrain,'*')

%% Poisson process with absolute refractory period -- from a random github page

% Real neurons are in a refractory state immediately after a spike. 
% Thus, during a certain ?refractory period? after a spike, it is 
% impossible or much harder to evoke the next spike.
% When the effect of refractoriness is taken into account, spikes 
% can no longer considered to be independent because the probability 
% of observing the next spike then also depends on the time that has 
% passed since the last spike.
%
% In point process models, an approximation of refractoriness is to 
% set the probability of generating a spike to 0 for a certain absolute 
% refractory period tr after a spike, and then to r afterwards. 
%
% We simulate a homogeneous Poisson process and add an absolute refractory
% period of t_r = 5 ms. We use different ?driving? rates 
% r = 10, 50, 100, 200, 500, and 1000Hz.
% We simulate 1000 spikes for each rate r.
%
% Adding an absolute refractory period corresponds to simply adding a 
% fixed value to each interspike interval of a ?normal? 
% homogeneous Poisson process.


rates_ref = [10 50 100 200 500 1000];  % "driving" rates in spikes per second
n = 1000;   % desired number of spikes
ISIs_ref = zeros(n,length(rates_ref));
randNumbers3 = rand(n,1);  % n random numbers, uniform distribution
t_r = 5;    % absolute refractory period in ms

for i = 1:length(rates_ref)
    ISIs_ref(:,i) = -log(randNumbers3)/rates_ref(i)*1000 + t_r;   % ISIs in milliseconds
end

% compute spike times from ISIs

SpikeTimes_ref = zeros(n,length(rates_ref));
for j = 1:length(rates_ref)
    SpikeTimes_ref(1,j) = ISIs_ref(1,j);
    for i = 2:n
        SpikeTimes_ref(i,j) = SpikeTimes_ref(i-1,j)+ISIs_ref(i,j);
    end
end