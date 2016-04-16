%% PCA on simulated neurons
%
% Guidance is provided via comments and example code below.
%
% To make this script easier to use, the task is broken into  sections,
% each section with a bold header (those starting with '%%') can be run by
% themselves (without running the entire script) using the 'run section
% command'. Just place your cursor within one of these sections (it will
% become highlighted) to allow this functionality.
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

% =======================
% Insert/Modify code here

psthFigName = 'psth_response_fig';
saveFormattedFig(psthFigName)

% =======================


%% Determine the number of neurons and number of time points 
% Use the 'size' function to get the size of each dimension of the data matrix 
%   hint: type 'doc size' for help on the size function
% Assign the number of neurons to the variable 'nNeurons' [ for n(umber)Neurons ]
% Assign the number of time points to the variable 'nTimePts'

% =======================
% Insert/Modify code here

% This is an 'empty' assignment, change it to be the appropriate values by
% using the size function.
nNeurons = []; 
nTimePts = [];

% =======================

%% "Center" the neural responses by subtracting off each mean response 
% First calculate the mean of each neuron's response using 'mean'. 
% Make sure you calculate the mean of data along the correct axis!
% Calculate the neuron's mean response and NOT the mean response of all
% neurons for one time point. 
%   hint: type 'doc mean' for how to choose the axis 'mean' averages over
%
% Assign the mean response matrix to the variable 'mu' below.

% =======================
% Insert/Modify code here

mu = [];

% =======================

% Next subtract each neuron's mean response from its response.
%   hint: The easiest way to do this is by subtracting the matrix mu from 
%       the data matrix. To perform subtraction the dimensions of the
%       matricies must match. You will need to make a copy of mu for
%       each time point to get the dimensions to be the same.
%   hint: 'repmat' will repeat a matrix to make a larger matrix, type 'doc
%       repmat' for help using this function. Use the nTimePts variable to
%       define the number of time dimensions you need
%
% Assign the mean subtracted (centered) responses to the variable
% 'centeredData' below.

% =======================
% Insert/Modify code here

centeredData = [];

% =======================

%% Plot centered neural responses for the first six neurons
% Use the plotting code above, or any method you want, to visualize mean
% subtraction. Include a title or otherwise indicate what operation has
% occured. You do not need to include the stimulus, but may if you like.

% =======================
% Insert/Modify code here

% =======================


%% Save a figure of the centered data
% Uncomment the code below and fill in a name for your figure.

% =======================
% Insert/Modify code here

%centeredFigName = '';
%saveFormattedFig(centeredFigName);

% =======================


%% Generate the covariance matrix
% Use the 'cov' function to calculate the covariance (using the centered
% data!). This will return a matrix that has variance along the diagonal
% entries and the covariance in the off-diagonal entries.

% You have done this correctly if dataCov(1,1) = XXXXXXXX

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


%% Find the principal components of the responses
% Use the 'pca' function on the data matrix (you do not need to use the
% centered data).
%   hint: Make sure the input to 'pca' is in the correct orientation, you
%       may need to transpose the data matrix using the transpose function 
%       e.g. x_transposed = x'
%
% The 'pca' function returns several values, for us the important ones are:
%   'score'     - essentially the "Principal Components" or PCs
%   'explained' - The amount of variance explained by each of the PCs
%   'coeff'     - also known as "loadings" these measure the importance of
%                 each variable in accounting for the variability of the 
%                 associated PC.
%                 (named coeff as they are the correlation coefficients 
%                 between the 'scores' and the original variables)
%
% Read the help documentation on pca for further information. Plotting
% values against each other can be helpful in exploring the data.

% You have done this correctly if coeff(1,1) = XXXXXXXX

% =======================
% Insert/Modify code here

% =======================

%% Plot explained variance (~Scree plot)
% Use the output from the 'pca' function above to make a plot of the
% different PC contributions to explained variance in the data.
%
%   hint: to make the plot more clear, pass '-o' to the 'plot' function.

% =======================
% Insert/Modify code here

% =======================

%% Save the explained variance figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================

%% Plot the first 6 principal components
% Use the output of 'pca' to plot the principal components. See the notes
% in the previous section for help on this.

% =======================
% Insert/Modify code here

% =======================

%% Save the Principal Component figure
% Use the example above to save the matrix.

% =======================
% Insert/Modify code here

% =======================

%% Find the covariance matrix of the principal components & plot it
% Remember this corresponds to the 'score' output of the 'pca' function.
% Again use the 'imagesc' and 'colorbar' functions for plotting.
% 
% You have done this correctly if the first entry in the matrix = XXXXXXXX

% =======================
% Insert/Modify code here

% =======================

%% Save the covariance matrix figure
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

%% Extension problems
% If you feel confident or would like to gain additional practice, please
% continue by answering the extension problems as outlined in the homework
% instructions. Some problems will require more coding please complete this
% below.


