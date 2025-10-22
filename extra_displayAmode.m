% this script is for showing the signal in the time domain without doing
% anything

clc; clear;

% path to our data
path_data        = 'data';
experiment       = 'ct1';
path_measurement = strcat(path_data, filesep, experiment);

path_outcomes    = 'outcomes';
path_labels      = strcat(path_outcomes, filesep, experiment);
filename_labels  = 'labeledSignalSet.mat';
fullpath_labels  = strcat(path_labels, filesep, filename_labels);

%% Load A-Mode Data

% add path for reading signal
path_usfunctions = 'functions\AModeMocap\';
addpath(path_usfunctions);

% get us data
dname = uigetdir(pwd);
[USData, ~, ~] = readTIFF_USsignal(dname, 30, 3500);

% preparing constants for data spesification
data_spec.n_ust     = size(USData, 1);
data_spec.n_samples = size(USData, 2);
data_spec.n_frames  = size(USData, 3);

% preparing constants for ultrasound spesification
us_spec.v_sound     = 1540; % m/s+
us_spec.sample_rate = 50 * 1e6; % Hz
us_spec.index2time_constant      = 1/(us_spec.sample_rate); %s (periode)
us_spec.index2distance_constant  = (1e3 * us_spec.v_sound) / (2 * us_spec.sample_rate);
us_spec.d_vector    = (1:data_spec.n_samples) .* us_spec.index2distance_constant; % mm
us_spec.t_vector    = ((1:data_spec.n_samples) .* us_spec.index2time_constant) * 1e6; % mu s

% signal processing to get m-mode data
path_barkercode = 'functions\AModeMocap\kenans_barkercode.txt';
[USsignals_envelop, USsignals_processed] = process_USsignal2(USData, data_spec, us_spec, path_barkercode, 'default');

% read the window
ust_config.WindowRange  = readmatrix(strcat(path_measurement, filesep, 'window_time.csv')) .* 1e-6; % in seconds
% convert windows in mm to windows in index
ust_config.WindowRange_i = floor(ust_config.WindowRange/us_spec.index2time_constant + 1);
% detect peaks
allpeaks    = peaks_USsignal_windowed3(USsignals_envelop, data_spec, us_spec, ust_config.WindowRange, ust_config.WindowRange_i, 'default');

% get the median of the peaks in respect to time
peak_time      = median(allpeaks.times, 2) * 1e6; % in microsecond
peak_width     = median(allpeaks.width, 2) * 1e6; % in microsecond
peak_sharpness = median(allpeaks.sharpness, 2) * 1e6; % in microsecond


%%

% prepare window
figure1 = figure('Name', 'A-mode Signal in Time Domain with Bone Surface Belief (manually annotated)', 'Position', [100 100 1800 700]);

% select frame
frame_selected = 25;
probe_selected = [20, 21, 24, 25, 26, 27, 28, 29, 30];

% variable to store uncertainties
uncertainty_time = zeros(data_spec.n_ust, 2);

% loop for selected frame
for frame_current = frame_selected
    
    % loop for selected probe
    for probe_current = probe_selected

        % get the current data
        USsignals_current = USData(probe_current, :, frame_current);
        USsignals_processed_current = USsignals_processed(probe_current, :, frame_current);
        USsignals_envelop_current   = USsignals_envelop(probe_current, :, frame_current);

        % make a figure with tabs
        tab = uitab('Title', num2str(probe_current));
        ax_amode  = axes(tab);

        % plot our ultrasound signal (raw and envelope) complete with all
        % necessary details for the plot (title, label, limit, etc).
        plot(ax_amode, us_spec.t_vector, USsignals_processed_current, '-g');
        hold(ax_amode, 'on');
        plot(ax_amode, us_spec.t_vector, USsignals_current, '-b');
        plot(ax_amode, us_spec.t_vector, USsignals_envelop_current, '-r', 'LineWidth', 3);
        title(ax_amode, 'A-mode Raw Signal', 'Interpreter', 'latex');
        axis(ax_amode, 'tight');
        xlabel(ax_amode, 'Time ($\mu$s)', 'Interpreter', 'Latex');
        ylabel(ax_amode, 'Amplitude', 'Interpreter', 'Latex');
        ax_amode.XGrid = 'on';
        ax_amode.XMinorGrid = 'on';
        ylim(ax_amode, [-12000, 12000]);

        % display the peak from the signal
        current_peaktime      = allpeaks.times(probe_current, frame_selected) * 1e6;
        current_peaksharpness = allpeaks.sharpness(probe_current, frame_selected);
        plot(ax_amode, current_peaktime, current_peaksharpness, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 9);

%         % display the peak width
%         current_peakwidth     = allpeaks.width(probe_current, frame_selected) * 1e6;
%         peakwidth_border      = current_peaktime + [-1 1] * (current_peakwidth/2);
%         xline(peakwidth_border(1), '-b');
%         xline(peakwidth_border(2), '-b', {'Peak Width'});

        % store the uncertainties
%         uncertainty_time(probe_current, :)      = [peakwidth_border(1) current_peaktime];
    
    end
end

%%
% 
% % fat and muscle
% uncertainty_vsound = [1440.2 1588.4];
% density            = 10;
% discreteDistribution_uncertaintyVsound       = linspace(uncertainty_vsound(1), uncertainty_vsound(2), density);
% discreteDistribution_uncertaintyDistribution = {}; 
% 
% figure2 = figure;
% 
% % loop for selected probe
% for probe_current = probe_selected
%     clf;
% 
%     discreteDistribution_uncertantyTime = linspace(uncertainty_time(probe_current,1), uncertainty_time(probe_current,2), density);
% 
%     uncertantyMatrix    = (discreteDistribution_uncertantyTime' * 1e-6) * (discreteDistribution_uncertaintyVsound*1e3);
%     uncertantyMatrix    = reshape(uncertantyMatrix, 1, density*density)';
% 
%     discreteDistribution_uncertaintyDistribution{probe_current} = fitdist(uncertantyMatrix, 'Normal');
% 
%     x_values = linspace(0, 60, 1000);
%     y_values = pdf(discreteDistribution_uncertaintyDistribution{probe_current}, x_values);
%     plot(x_values, y_values);
%     hold on;
%     plot(uncertantyMatrix, zeros(size(uncertantyMatrix)), 'x');
%     xlabel('Depth (mm)');
%     ylabel('Likelihood');
% 
% 
% end

% test = reshape(discreteDistribution_uncertaintyDistance(25, :, :), 1, density*density)';
% pd   = fitdist(test, 'Normal');
% 
% x_values = linspace(0, 20, 1000);
% y_values = pdf(pd, x_values);
% figure;
% plot(x_values, y_values);
% hold on;
% plot(test, zeros(size(test)), 'x');































