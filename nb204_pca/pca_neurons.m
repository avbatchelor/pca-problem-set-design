%% PCA on simulated neurons
%
% Guidance is provided via comments and example code below.
%
% To make this script easier to use, the task is broken into  sections,
% each section with a bold header (those starting with '%%') can be run by
% themselves (without running the entire script) using "ctrl + enter"
% (windows) or "command + enter" (MAC). Just place your cursor within one of
% these sections (it will become highlighted) to allow this functionality.
% 
% AVB & SLH 4/2016

%% Close figures and clear workspace
clear;      % Delete all  variables in workspace, you will lose unsaved variables
close all;  % Close all of the open figure windows

%% Load data
% This line loads three variables: data, stim, time
load('pca_data.mat')

% data is a 58x5000 matrix, Neurons x Time Points
% Each row is the PSTH of a neuron's response to the stimuli
% stim is a 1x5000 column vector of the Stimuli over time
% time is a 1x5000 column vector of the Time in second

%% Plot data for six neurons
% This section of the code plots the stimulus and the responses for the
% first six neurons in the same figure.  You can copy and change this code
% to create your other figures.

% If you can't see the data in the figure, maximize the figure so you can. 

figure
ax(1) = subplot(6,1,1); % subplot allows you to plot multiple graphs in the same figure
plot(time,stim','r') % Plot the stimulus in red ('r')
ylabel('Odor concentration')
title('Stimulus')
for i = 1:5 % Loop over first 5 neurons to plot their responses
    ax(i+1) = subplot(6,1,i+1);
    plot(time,data(i,:)) % This plots row i of data
    title(['Response of neuron ',num2str(i)])
end
xlabel('Time (seconds)')
ylabel('Spike rate (Hz)')

%% Save figure 
% The saveFormattedFig function saves the last figure you have made in a
% format that is easy to insert (simply by dragging) into a word document.
% You can use this code to save all your other figures by calling
% 'saveFormattedFig' with the the name of the file you wish to save as a
% string, as below:

% You don't need to change any code in the saveFormattedFig.m file. 
psthFigName = 'psth_response_fig';
saveFormattedFig(psthFigName)

%% Generate the covariance matrix
% Use the 'cov' function to calculate the covariance (using the centered
% data!). This will return a matrix that has variance along the diagonal
% entries and the covariance in the off-diagonal entries.

% You have done this correctly if dataCov(1,1) = 0.3275
% Replace the [] with your own code below.

% =======================
% Insert/Modify code here

dataCov = [];

% =======================

%% Plot the covariance matrix
% You can use the 'imagesc' function to visualise the covariance matrix.
% Calling the 'colorbar' function adds the color scale on your graph.
% 

% =======================
% Insert/Modify code here

% =======================

%% Save a covariance matrix figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================

%% Perform PCA
% Use the 'pca' function to run PCA on the data matrix. Your goal is to
% analyze the relationships among neurons, not the relationships among
% time points. Your hypothesis is that the responses of all 58 neurons
% can be reduced to linear combinations of a few orthogonal basis
% functions, where each basis function can be conceptualized as a
% distinct "response type". You should set up the PCA so that you obtain
% 58 PCs; you hypothesize that only a few of these PCs are needed to
% explain most of the variance in the data.

% The 'pca' function returns several values. For us the important ones are:
% 
%   'score' - This is a matrix containing the representation of the data
%       set in PC space. Essentially, this is the data set after it has
%       been rotated. Recall that the original data set consisted of 58 
%       vectors, with each vector representing a neural response measured 
%       at 5000 time points; this new matrix therefore also consists of 58 
%       vectors measured at 5000 time points. The vectors in the score 
%       matrix are also sometimes referred to as "PCs".
%   'explained' - This is a list of numbers quantifying the percentage of
%       the variance in the data explained by each of the PCs, in 
%       descending order of variance explained.
%   'coeff' - This is a matrix that quantifies the importance of each
%       variable (here, each neuron in the original dataset) in accounting 
%       for the variability of the associated PC. (This matrix gets the 
%       name 'coeff' because it contains the correlation coefficients 
%       between the data in PC space and the original data.) These values 
%       are also known as loadings. Each column of coeff contains 
%       coefficients for one principal component, and the columns are in 
%       descending order of component variance.

% Read the help documentation on pca for further information. 

% If you're confused, we recommend first performing the plotting steps
% below. Plotting the outputs from the pca function can be helpful for
% understanding what pca is doing.

% You have done this correctly if coeff(1,1) = 0.0724

% =======================
% Insert/Modify code here

% =======================

%% Use the output of 'pca' to plot the first six principal components 
% The first six PCs are the first six vectors in the score matrix, where
% each vector is a list of 5000 numbers.

% =======================
% Insert/Modify code here

% =======================

%% Save the Principal Component figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================

%% Plot explained variance (~Scree plot)
% Use the output from the 'pca' function above to make a plot of the
% different PC contributions to explained variance in the data.
%
%   hint: To make it easy to see the variance explained by each pc when you
%       plot 'explained' also pass '-o' to the plot function, like this
%       example: plot(explained,'-o')

% =======================
% Insert/Modify code here

% =======================

%% Save the explained variance figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================

%% Make a 3D plot of each neuron's loadings for the first 3 PCs
% Use 'plot3' to make a 3D plot. Plot each loading as a discrete dot or
% circle for clarity, and please label the axes. Remember the loadings
% correspond to the 'coeff' output of the 'pca' function.
%   hint: type 'doc Chart Line Properties' and select the first search
%       result for help with plot formatting.
%   hint: using 'grid on' might make your graph more easily viewable

% =======================
% Insert/Modify code here

% =======================

%% Save loadings figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================

%% Find the covariance matrix of data in the PC space and plot it
% Use the 'imagesc' function and the 'colorbar' function for plotting. You
% should have 58 PCs, so this should be a 58-by-58 matrix (i.e. the same
% size as the previous covariance matrix you plotted).
% 
% You have done this correctly if the first entry in the matrix = 47.1669

% =======================
% Insert/Modify code here

% =======================

%% Save the covariance matrix figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================


%% Extension problems
% If you feel confident or would like to gain additional practice, please
% continue by answering the extension problems as outlined in the homework
% instructions. Some problems will require more coding, you may complete
% this below.

