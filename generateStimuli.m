function [stimuli,time] = generateStimuli(sampleRate)
% stimuli [NxM Matrix] N = stimulus type, M = time in ms
% stimuli are lowpass filtered and rectified to make them more "realistic"
% stimuli generated/returned at sampleRate
%
% time - timebase returned in seconds
%
% SLH 2016

F = sampleRate/1000;

% Lowpass filter for more realistic stimuli

% Analog 2kHz 5th order bessel filter
[z,p,k] = besself(5,2000);
% Analog to digital mapping
[zd,pd,kd ] = bilinear(z,p,k,sampleRate);
% Convert to SOS form
[sos] = zp2sos(zd,pd,kd);

% matrix of stim[ms]
stimuli = zeros(5,1000*F);

for stimType = 1:6
    switch stimType
        case 1
            % 500ms on
            stimCmd = [0.*ones(250*F,1); ones(500*F,1); 0.*ones(250*F,1)];
        case 2
            % 100ms on/off
            cycle = [repmat([ones(100*F,1); 0.*ones(100*F,1)],2,1); ones(100*F,1)];
            stimCmd = [0.*ones(250*F,1); cycle ; 0.*ones(250*F,1)];
        case 3
            % 50ms on/off
            cycle = repmat([ones(50*F,1); 0.*ones(50*F,1)],5,1);
            stimCmd = [0.*ones(250*F,1); cycle ; 0.*ones(250*F,1)];
        case 4
            % ramping up over 500ms
            ramp = linspace(0,1,500*F);
            stimCmd = [0.*ones(250*F,1); ramp'; 0.*ones(250*F,1)];
        case 5
            % ramping down over 500ms
            ramp = linspace(1,0,500*F);
            stimCmd = [0.*ones(250*F,1); ramp'; 0.*ones(250*F,1)];
        case 6
            % blank
            stimCmd = 0.*ones(1000*F,1);
    end

    % filter, rectify, normalize stimuli
    stimuli(stimType,:) = sosfilt(sos,stimCmd);
    stimuli(stimType,[stimuli(stimType,:) < 0]) = 0;
    if max(stimuli(stimType,:)) > 0
        stimuli(stimType,:) = stimuli(stimType,:)./max(stimuli(stimType,:));
    end
    time = (1:1000*F)./sampleRate;
end