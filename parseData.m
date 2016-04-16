%% Load unparsed data 
load('C:\Users\Alex\Documents\GitHub\nb204-LNLNP-tuning-simulation\blind_neurons.mat')

%% Make data matrix 

% data are already shuffled? why is this like this
allN_cell = struct2cell(allN);
allN_mat = cell2mat(allN_cell);
allN_mat_p = permute(allN_mat,[2 1 3]);
N_mat = reshape(allN_mat_p,6000,58);
N_mat = N_mat';

data = N_mat;

all_stim = reshape(stim',[],1);
stim_down = 5*downsample(all_stim,5);
samprate = 1/1000;
time = samprate:samprate:6;
stim = stim_down; 
stim = stim';

save('C:\Users\Alex\Documents\GitHub\nb204-LNLNP-tuning-simulation\pca_data.mat','data','stim','time')
