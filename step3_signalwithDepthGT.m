%% SUMMARY
% This code is to generate the depth estimated by the ultrasound signal 
% (from user selection) and the ground truth.
%
% Ultrasound signal was stored and ecoded in an .TIFF file. Find the whole
% data in: 'data_organized\US\CTexp1\duringCT_cut' or
%          'data_organized\US\CTexp2\duringCT_cut'
% To read this data use the function readTIFF_USsignal(). You can find the
% function in: functions\AModeMocap\readTIFF_USsignal.m
%
% User selection is defined in the experiment and stored in an *.ini file.
% Find it in: 'data_organized\US\CTexp1' or 
%             'data_organized\US\CTexp2'

% IMPORTANT NOTE: 
% Please specify the correct folder name for path_measurement. You can 
% choose 'ct1' or 'ct2'.

clc; clear;

% path to our data
path_data        = 'data_organized';
% path_measurement = strcat(path_data, filesep, 'CT', filesep, 'ct1');
path_measurement = strcat(path_data, filesep, 'CT', filesep, 'ct2');

% load ground truth data (tube based)
filename_transducermat = 'transducers2_probebased_tubedepth.mat';
fullpath_transducermat = strcat(path_measurement, filesep, filename_transducermat);
load(fullpath_transducermat);
transducers_tubedepth = transducers;

% load ground truth data (line based)
filename_transducermat = 'transducers2_probebased_linedepth.mat';
fullpath_transducermat = strcat(path_measurement, filesep, filename_transducermat);
load(fullpath_transducermat);
transducers_linedepth = transducers;

%% A-mode signal

% add path for reading signal
path_usfunctions = 'functions\AModeMocap\';
addpath(path_usfunctions);

% get us data
fprintf('Select the Ultrasound Signal folder\n');
dname = uigetdir(pwd);
[USData, ~, ~] = readTIFF_USsignal(dname, 30, 3500);

% preparing constants for data spesification
data_spec.n_ust     = size(USData, 1);
data_spec.n_samples = size(USData, 2);
data_spec.n_frames  = size(USData, 3);

% preparing constants for ultrasound spesification
us_spec.v_sound     = 1540; % m/s
us_spec.sample_rate = 50 * 1e6; % Hz
us_spec.period      = 1/(us_spec.sample_rate); %s
us_spec.index2distance_constant  = (1e3 * us_spec.v_sound) / (2 * us_spec.sample_rate);
us_spec.d_vector    = (1:data_spec.n_samples) .* us_spec.index2distance_constant; % mm
us_spec.t_vector    = ((1:data_spec.n_samples) .* us_spec.period) * 1e6; % mu s

% Signal processing to get envelopes data. Note: i give us_spec as the
% argument for this function, even thought us_spec has a lot of fields
% (initialized above) but in this function i only used us_spec.sample_rate
path_barkercode = 'functions\AModeMocap\kenans_barkercode.txt';
[USsignals_envelop, USsignals_processed] = process_USsignal2(USData, data_spec, us_spec, path_barkercode, 'default');


%% Main Loop

% prepare window
figure1 = figure('Name', 'Signal with GT (Without Velocity Correction)', 'Position', [100 100 1800 600]);

% select frame, here i only use one frame, because i want to display. if
% you want to use all frame (let's say you want to see the mean/median of
% the results), you should disable the display.
frame_selected = 10;
probe_selected = [20 21 24 25 26 27 28 29 30];

% change to false if you want to show time (mu s) instead of depth (mm)
isdepth = true;

% loop for selected frame
for frame_current = frame_selected
    
    % loop for selected probe
    for probe_current = probe_selected

        % get the current ultrasound signal
        USsignals_processed_current = USsignals_processed(probe_current, :, frame_current);
        USsignals_envelop_current   = USsignals_envelop(probe_current, :, frame_current);

        % make a figure with tabs
        tab = uitab('Title', num2str(probe_current));
        ax_amode  = axes(tab);

        % plot our ultrasound signal (raw and envelope) complete with all
        % necessary details for the plot (title, label, limit, etc).
        if(isdepth)
            plot(ax_amode, us_spec.d_vector, USsignals_processed_current, '-g');
            hold(ax_amode, 'on');
            plot(ax_amode, us_spec.d_vector, USsignals_envelop_current, '-r', 'LineWidth', 2);
        else
            plot(ax_amode, us_spec.t_vector, USsignals_processed_current, '-g');
            hold(ax_amode, 'on');
            plot(ax_amode, us_spec.t_vector, USsignals_envelop_current, '-r', 'LineWidth', 2);
        end
        axis(ax_amode, 'tight');
        title(ax_amode, 'A-mode Raw Signal', 'Interpreter', 'latex');
        if(isdepth)
            xlabel(ax_amode, 'Depth (mm)', 'Interpreter', 'Latex');
        else
            xlabel(ax_amode, 'Time ($\mu$s)', 'Interpreter', 'Latex');
        end
        ylabel(ax_amode, 'Amplitude', 'Interpreter', 'Latex');
        ax_amode.XGrid = 'on';
        ax_amode.XMinorGrid = 'on';
        ylim(ax_amode, [-12000, 12000]);

        % We want to plot the depth, but some of the probes can detect
        % nothing. To avoid the error, let's check our data first, if the
        % depth is empty or not. Here i use depth.max, but other value,
        % such as depth.min or depth mean can also be used.
        if(~isempty(transducers(probe_current).depth.max) && isdepth)

            % show the min and the max by using line
            xline(ax_amode, transducers_tubedepth(probe_current).depth.min, '-b', 'Tag', 'plot_amode', 'LineWidth', 1);
            xline(ax_amode, transducers_tubedepth(probe_current).depth.max, '-b', 'Tag', 'plot_amode', 'LineWidth', 1);

            % this if block below is just for display purposes, to let the
            % matlab know where should i put the label
            if(transducers_tubedepth(probe_current).depth.mean < transducers_linedepth(probe_current).depth.mean)
                tubedepth_textposition = 'left';
                linedepth_textposition = 'right';
            else
                tubedepth_textposition = 'right';
                linedepth_textposition = 'left';
            end

            % show the depth from tube method
            xline( ax_amode, transducers_tubedepth(probe_current).depth.mean, ...
                   '--m', {num2str(transducers_tubedepth(probe_current).depth.mean)}, ...
                   'LabelHorizontalAlignment', tubedepth_textposition, ...
                   'Tag', 'plot_amode', 'LineWidth', 1);

            % show the line from line method
            xline( ax_amode, transducers_linedepth(probe_current).depth.mean, ...
                   '-.', {num2str(transducers_linedepth(probe_current).depth.mean)}, ...
                   'LabelHorizontalAlignment', linedepth_textposition, ...
                   'Color', '#7E2F8E', 'Tag', 'plot_amode', 'LineWidth', 1);
        end

        xlim(ax_amode, [0 35]);

        % display the legend
        legend({'A-mode Signal', 'Envelope', 'Min Distance', 'Max Distance', 'Mean Distance', 'Line Intersection'});
    
    % end probe loop 
    end

% end frame loop
end




