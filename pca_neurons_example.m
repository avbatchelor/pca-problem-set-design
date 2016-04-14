%% PCA on simulated neurons

% AVB 3/10/2016

close all 
clear all

%% Load data 
load('C:\Users\Alex\Documents\GitHub\nb204-LNLNP-tuning-simulation\pca_data.mat')

%% Plot data
figure
ax(1) = subplot(7,1,1);
plot(time,stim','r')
ylabel('Odor concentration')
title('Stimulus')
for i = 1:6
    ax(i+1) = subplot(7,1,i+1);
    plot(time,data(i,:))
    title(['Response of neuron ',num2str(i)])
end
xlabel('Time (seconds)')
ylabel('Spike rate (Hz)')
linkaxes(ax(:),'xy')
ax_lims = [min(data(:)),max(data(:))];
ylim(ax_lims)

%% Get the number of neurons and time points
[num_neurons,num_time_points] = size(data);

%% Subtract off the mean of each feature 
mu = mean(data,2);
data_cent = data - repmat(mu,1,num_time_points);
new_mu = mean(data_cent,2);


%% Plot centered data
figure
plot(time,data_cent)
title('centered data')

%% Calculate the covariance matrix and plot
cov_mat = (1/num_neurons).*data_cent*data_cent';
% Plot covariance matrix 
figure
imagesc(cov_mat)
colorbar
title('covariance matrix')

%% Find the eigenvectors and eigenvalues 
[V,D] = eig(cov_mat); % columns of V are eigenvectors 
% Make a vector of eigenvalues 
lambda = diag(D);
% Sort eigenvalues with largest first 
[lambda_sorted, Idx] = sort(lambda,'descend');

% plot the eigenvalues 
figure
plot(lambda_sorted,'o')
title('eigenvalues')

%% Plot the eigenvectors
figure
title('eigenvectors')
for i = 1:6
    ax(i) = subplot(6,1,i);
    plot(V(:,Idx(i)))
end
linkaxes(ax(:),'y')
ax_lims = [min(V(:)),max(V(:))];
ylim(ax_lims)

%% Find principal components by projecting data onto the eigenvectors 
pcs = data'*V;
% Plot principal components
figure
title('principal components')
for i = 1:6
    ax(i) = subplot(6,1,i);
    plot(pcs(:,i))
end
linkaxes(ax(:),'y')
ax_lims = [min(pcs(:)),max(pcs(:))];
ylim(ax_lims)


%% Find weights by projecting each cell onto principal components
% Find weights
wts = data_cent*pcs;
% Plot weights
figure
plot3(wts(:,1),wts(:,2),wts(:,3),'o')
xlabel('projection onto PC1')
ylabel('projection onto PC2')
zlabel('projection onto PC3')
title('Projection of neurons onto PCs 1, 2 and 3')