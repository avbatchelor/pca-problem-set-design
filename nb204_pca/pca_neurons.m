%% PCA on simulated neurons

% AVB & SLH 3/10/2016

%% Close figures and clear workspace
close all 
clear all 
% Note that this command clears the workspace so you will lose any
% variables you haven't saved

%% Set some defaults for figure formatting so your figures look nice 
set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultFigureColor','w')
figure 
close all 

%% Load data 
load('pca_data.mat')
% This line loads three variables: data, stim, time. 
% data is a 58 x 6000 matrix because there are 58 neurons and 6000 time points.  
% Each row is the PSTH of a simulated neuron in response to the stimuli. 
% stim is a 1 x 6000 vector of the stimuli over time. 
% time is a 1 x 6000 vector of the time in seconds.

%% Plot data for the first six neurons
% This section of the code plots the stimulus and the responses for the
% first six neurons in the same figure.  You can copy and change this code
% to create your other figures. 
figure
ax(1) = subplot(7,1,1); % subplot allows you to plot multiple graphs in the same figure
plot(time,stim','r') % Plot the stimulus
ylabel('Odor concentration')
title('Stimulus')
for i = 1:6 % Loop through the first six neurons and plot their responses
    ax(i+1) = subplot(7,1,i+1);
    plot(time,data(i,:)) % This plots row i of data
    title(['Response of neuron ',num2str(i)])
end
xlabel('Time (seconds)')
ylabel('Spike rate (Hz)')
linkaxes(ax(:),'xy') % This makes sure the axes have the same scaling 
ax_lims = [min(data(:)),max(data(:))]; % Finds the correct axis limits so you can see all the data 
ylim(ax_lims) % Set the axis limits 

%% Save the figure 
% This code saves the latest figure you have made in a format that is easy
% to insert into a word document.  You can use this code to save all your
% other figures by inserting it after the code that makes the figure. You
% just need to specify the name of the file as a string.  In this case the
% name is 'neuron_data_fig'.
my_save_fig('neuron_data_fig')

%% Determine the number of neurons and number of time points 
% You can use the function 'size' to find the size of each dimension of the
% data matrix 

% ================
% Insert code here
% ================

%% Subtract off the mean response of each neuron (This is 'centering')
% First calculate the mean of each neuron's response using the 'mean'
% function.  Make sure you calculate the mean of data along the correct axis
% (i.e. calculate the neuron's mean response and not the mean response of
% all neurons for one time point). You can specify the axis as an input
% argument for the mean function.

% ================
% Insert code here
% ================

% Next subtract the mean response from the data.  Note when you subtract
% two matrices they need to have the same dimensions.  You can use 'repmat'
% to make a matrix of mean values. 

% ================
% Insert code here
% ================

%% Plot centered data for the first six neurons

% ================
% Insert code here
% ================


%% Save a figure of the centered data

% ================
% Insert code here
% ================


%% Generate the covariance matrix
% Please do not use MATLAB's 'cov' function for generating the covariance
% matrix but instead generate it by multiplying the appropriate matrices as
% we discussed in class. You can multiply to matrices A and B using
% mtimes(A,B) or A*B

% Think about what the dimensions of your covariance matrix should be. If
% your original data has M rows and N columns and you are performing
% dimensionality reduction to reduce the number of rows then your covarianc
% matrix should be an M x M matrix.  

% To check that you've calculated the covariance matrix correctly, the
% first three values of the first row should be 327.7398  123.2678
% 302.6009.

% ================
% Insert code here
% ================


%% Plot the covariance matrix
% You can use the imagesc function to visualise the covariance matrix. 
% The colorbar function puts the color scale on your graph. 

% ================
% Insert code here
% ================

%% Save the covariance matrix figure

% ================
% Insert code here
% ================


%% Find the eigenvectors and eigenvalues 
% Use the 'eig' function to do this. 

% ================
% Insert code here
% ================

%% Make a vector of eigenvalues
% Use the 'diag' function to extract the eigenvalues from the matrix output
% of eig. 

% ================
% Insert code here
% ================

%% Sort the eigenvalues and eigenvectors 
% Use the 'sort' function to sort the eigenvalues from highest to lowest 
% Make sure to get the indexes as an output from the sort function so that
% you can sort the eigenvectors too. 

% ================
% Insert code here
% ================

% Now sort the eigenvectors so that they are in the same order as the
% eigenvalues.  To do this don't use the 'sort' function but instead index
% the eigenvector matrix using the indexes you got from sorting the
% eigenvalues. 

% ================
% Insert code here
% ================


%% Plot the eigenvalues 
% Plot the eigenvalue value vs. the eigenvalue number. Plot each eigenvalue
% as a discrete dot rather than a line that joins each eigenvalue so that
% it is easier to see the value for each eigenvalue. 

% ================
% Insert code here
% ================

%% Save eigenvalues plot

% ================
% Insert code here
% ================

%% Plot the first six eigenvectors

% ================
% Insert code here
% ================

%% Save eigenvectors plot

% ================
% Insert code here
% ================

%% Find principal components by projecting data onto the eigenvectors 
% Recall that you project data onto the eigenvectors by multiplying the
% data matrix and the eigenvector matrix together.

% Again think about what dimensions your principal components should have.
% You may need to transpose your data matrix using the 'transpose'
% function. 

% ================
% Insert code here
% ================

%% Plot the first six principal components

% ================
% Insert code here
% ================

%% Save principal components plot

% ================
% Insert code here
% ================

%% Find the loadings of each neuron on the principal components
% You find these loadings by projecting the data matrix onto the principal
% components.  You do this my multiplying the data matrix and the principal
% components matrix.  Again think about the dimensions that the loadings
% should have. 

% ================
% Insert code here
% ================

%% Make a 3D plots of the loadings of each neuron on the first three principal components
% Use 'plot3' to make a 3D plot
% Make sure to plot each loading as a discrete dot. 
% Remember to label the axes

% ================
% Insert code here
% ================

%% Save loadings plot 

% ================
% Insert code here
% ================

