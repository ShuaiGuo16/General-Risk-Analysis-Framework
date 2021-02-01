clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Obtain FDF phase uncertainty description at a 26-by-10 
%            frequency-amplitude grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%   ===> (1) Select only limited phase measurements from full experimental 
%            results, assume 5% uncertainty for phase measurements
%   ===> (2) Use stochastic Gaussian Process regression to interpolate 
%            the selected uncertain measurements to a finer grid (26-by-10)
%   ===> (3) Record the mean and covariance of the phase data at each 
%            location of the finer grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: Optimization, Kriging scripts provided in the companion code of [1]
% Ref: [1] A. Forrester, Engineering Design via Surrogate Modelling: A Practical Guide
%             2008, Wiley.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select limited FDF phase data from experimental results
load './data/FDF_A_ori.mat'       % Full experimental data
Amp = [0.07;0.15;0.3;0.41;0.51;0.71];
Freq_list = [0;30;60;80;140;170;190;210;230;250];   
sliced_phase = Phase(ismember(Freq,Freq_list),:);
com_phase = Phase(~ismember(Freq,Freq_list),:);

%% Data visualization
figure(1)
% Training sample
Phase = sliced_phase'; 
[coor_freq,coor_amp] = meshgrid(Freq(ismember(Freq,Freq_list)),Amp);
scatter3(coor_freq(:),coor_amp(:),-Phase(:),'filled','MarkerFaceColor','r')
hold on
% Full sample
Phase = com_phase'; 
[coor_freq,coor_amp] = meshgrid(Freq(~ismember(Freq,Freq_list)),Amp);
scatter3(coor_freq(:),coor_amp(:),Phase(:),'filled','MarkerFaceColor','k')

%% Train the interpolation model
% Re-configure the sliced data
training_X = mesh2array(Freq(ismember(Freq,Freq_list)),Amp);
training_Y = sliced_phase';
training_Y = training_Y(:);

% Normalization
training_X(:,1) = training_X(:,1)/max(Freq_list);

% Set standard deviation value for each measurement 
phase_std = 1e-5*ones(size(training_Y,1),1);   % no uncertainty at 0Hz
phase_std(7:end) = 0.05*training_Y(7:end);       % uncertainty level: 5% of the mean

disp('Start training GP model')
[GP_Model] = GP_noise(training_X,training_Y,phase_std.^2);
disp('Training finished!')

%% Realization generation 
Freq_realization = 0:10:max(Freq_list);        % 26 levels of frequency
Amp_realization = linspace(0.07,0.71,10);    % 10 levels of amplitude
pred_X = mesh2array(Freq_realization',Amp_realization');
pred_X(:,1) = pred_X(:,1)/max(Freq_list);
[f,C] = pred_noise(pred_X, GP_Model);   % Mean and covariance matrix
% save './data/realization_phase.mat' GP_Model f C Freq_realization Amp_realization

%% Realization visualization
bootstrap = 5;
random_training_data = mvnrnd(f,C,bootstrap);
for i = 1:bootstrap
    disp('Start bootstrapping training')
    boot_f = random_training_data(i,:);   boot_f = boot_f';
    % Plot routine
    [coor_freq,coor_amp] = meshgrid(Freq_realization',Amp_realization');
    boot_f = reshape(boot_f,[length(Amp_realization),length(Freq_realization)]);
    surf(coor_freq,coor_amp,-boot_f,'EdgeColor','none','FaceColor','k','FaceAlpha',0.4)
end
hold off

xlabel('Frequency')
ylabel('Amplitude')
zlabel('Phase')
xticks(0:50:250)
yticks(0:0.2:0.8)
h = gca;
h.FontSize = 14;
view(16,28)