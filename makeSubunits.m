function [filts, weights] = makeSubunits(props,sampleRate)
% return 6 filters == "subunits" with weights from the props struct

if props.use(1)
    %%%% Filter 1: positive of differential ("ON") + small integration
    fltLen = (25/1000)*sampleRate;          % Sample width samples in filter
    hwDg = 1;                               % Half-width of filter
    tFilt = linspace(-hwDg,hwDg,fltLen);    % filter timebase
    stdGaus = 0.8;                          % std of Gaussian, in deg.
    gauss = exp(-tFilt.^2/stdGaus.^2);      % Gaussian envelope
    
    b = 4/1000;
    flt = ( +5 * tFilt + b ) .* gauss;
    flt = flt./max(flt);
    
    % Get weight from distribution
    weights{1} = (props.sigma(1) * rand(1)) + props.mu(1);
    filts{1} = flt;
else
    weights{1} = 0;
    filts{1} = 1;
end

%%%% Filter 2: negative of differential ("OFF") + small integration
if props.use(2)
    fltLen = (25/1000)*sampleRate;          % Sample width samples in filter
    hwDg = 1;                               % Half-width of filter
    tFilt = linspace(-hwDg,hwDg,fltLen);    % filter timebase
    stdGaus = 0.8;                          % std of Gaussian, in deg.
    gauss = exp(-tFilt.^2/stdGaus.^2);      % Gaussian envelope
    
    b = 4/1000;
    flt = ( -5 * tFilt + b ) .* gauss;
    flt = flt./max(flt);

    % Get weight from distribution
    weights{2} = (props.sigma(2) * rand(1)) + props.mu(2);
    filts{2} = flt;
else
    weights{2} = 0;
    filts{2} = 1;
end

%%%% Filter 3: integrating filter F = (alpha +/- raised cosines)
if props.use(3)

    fltLen = (50/1000)*sampleRate;
    tFilt = linspace(0,1,fltLen);
    tau = 70/1000;
    alphaF = tFilt.*exp(-tFilt/tau);
    alphaF = 5*alphaF./max(alphaF);
    cosF1 = (cos(linspace(-pi/2-(pi/12),pi/2-(pi/12),fltLen)));
    cosF2 = (cos(linspace(-pi/2,pi/2,fltLen)));
    combF = alphaF-.5*cosF1+2.5*cosF2;
    combF = combF./max(combF);

    % Get weight from distribution
    weights{3} = (props.sigma(3) * rand(1)) + props.mu(3);
    filts{3} = combF;
else
    weights{3} = 0;
    filts{3} = 1;
end

%%%% Filter 4: funky differentiating filter F = (alpha +/- raised cosines)
if props.use(4)
    %%
    fltLen = (45/1000)*sampleRate;
    tFilt = linspace(0,1,fltLen);
    tau = 20/1000;
    alphaF = tFilt.*exp(-tFilt/tau);
    alphaF = -alphaF./max(alphaF);
    cosF1 = (sin(linspace(-pi/2-(pi/6),pi/2-(pi/6),fltLen)));
    cosF2 = (cos(linspace(-pi/2+(pi/6),pi/2+(pi/6),fltLen)));
    b = 3000/1000;
    combF = 30*alphaF-3*cosF1+6*cosF2 - b;
    combF = combF./max(abs(combF));
    %%
    % Get weight from distribution
    weights{4} = (props.sigma(4) * rand(1)) + props.mu(4);
    filts{4} = combF;
else
    weights{4} = 0;
    filts{4} = 1;
end

%%%% Filter 5: triphasic slightly integrating filter
if props.use(5)
    %%
    fltLen = ((15)/1000)*sampleRate;   % desired filter length
    fltLen2 = ((31)/1000)*sampleRate;  % upsampled filter
    hwDg = 2;                          % Half-width of filter
    tFilt = linspace(-hwDg,hwDg,fltLen2); % upsampled filter timebase
    stdGaus = 0.75;                       % std of Gaussian, in deg.
    gauss = exp(-tFilt.^2/stdGaus.^2);   % Gaussian envelope
    mHatF = diff(diff(gauss));         % make upsampled mexican hat
    mHatF = downsample(mHatF,2);       % downsample mexican hat
    mHatF = mHatF(1:fltLen);           % truncate to desired length
    mHatF = mHatF./-(min(mHatF));
    expF = (1/50).^linspace(0,1,fltLen);
    expF = expF./max(expF);
    combF = -(2*expF .* mHatF - mHatF);
    combF = combF./max(combF);
    
    %%
    % Get weight from distribution
    weights{5} = (props.sigma(5) * rand(1)) + props.mu(5);
    filts{5} = combF;
else
    weights{5} = 0;
    filts{5} = 1;
end

%%%% Filter 6: triphasic slightly integrating filter
if props.use(6)
    %%
    fltLen = ((40)/1000)*sampleRate;   % desired filter length
    fltLen2 = ((81)/1000)*sampleRate;  % upsampled filter
    hwDg = 2;                          % Half-width of filter
    tFilt = linspace(-hwDg,hwDg,fltLen2); % upsampled filter timebase
    stdGaus = 0.75;                       % std of Gaussian, in deg.
    gauss = exp(-tFilt.^2/stdGaus.^2);   % Gaussian envelope
    mHatF = diff(diff(gauss));         % make upsampled mexican hat
    mHatF = downsample(mHatF,2);       % downsample mexican hat
    mHatF = mHatF(1:fltLen);           % truncate to desired length
    mHatF = mHatF./-(min(mHatF));
    expF = (20).^linspace(0,1,fltLen);
    expF = expF./max(expF);
    m = -20/1000;
    b = -110/1000;
    combF = +(4.*expF .* 2.*mHatF) + (m.*tFilt(1:fltLen) + b);
    combF = combF./max(combF);

    % Get weight from distribution
    weights{6} = (props.sigma(6) * rand(1)) + props.mu(6);
    filts{6} = combF;
else
    weights{6} = 0;
    filts{6} = 1;
end