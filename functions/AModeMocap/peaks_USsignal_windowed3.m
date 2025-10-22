function allpeaks = peaks_USsignal_windowed3(envelope_data, data_spec, us_spec, windowrange, windowrange_i, processmode)
% SUMMARY: The different between this function and other version is how do 
% we clip the data. peaks_USsignal_windowed - peaks_USsignal_windowed2 clip 
% the data by clippping the raw signal, then process the raw signal. Here, 
% we clip after we process the raw signal. The reason is, when the window is
% to small, some signal processing (filtering/correlation) can heavily 
% change the information of the raw data.
% The input of this function is the processed data, that is the envelope.

% to see the time
t_dsp = zeros(data_spec.n_frames, 1);

% put indicator to terminal
disp("Peak Detection is running, please wait ...");

% show the progress bar, so that the user is not bored
progress_bar = waitbar(0, sprintf('%d/%d Frame', 0, data_spec.n_frames), 'Name', 'Running Peak Detection');

% remove the warning if the findpeak function can't find the peak below
% threshold
warning('off', 'signal:findpeaks:largeMinPeakHeight');

for j=1:data_spec.n_frames
    
    % display progress bar
    if (mod(j,25)==0)
        waitbar( j/data_spec.n_frames, progress_bar, sprintf('%d/%d Frame', j, data_spec.n_frames) );
    end
    
    tic;
    for i=1:data_spec.n_ust

        % clip the envelope data
        data_clipped    = envelope_data(i, windowrange_i(i,1):windowrange_i(i,2), j);

        % different parameter for different envelope
        if(strcmp(processmode, 'default'))
            minpeakheight     = 100;
            minpeakprominence = 300;
        elseif(strcmp(processmode, 'cwt1') || strcmp(processmode, 'cwt2') || strcmp(processmode, 'cwt'))
            minpeakheight     = 100;
            minpeakprominence = 130;
        end

        % find local maxima
        [peaks, locs, width, prominence] =  findpeaks( data_clipped, 'MinPeakHeight', minpeakheight, 'MinPeakProminence', minpeakprominence, 'SortStr', 'descend', 'WidthReference', 'halfprom');

        % we only store the locs value if it is not empty, or else, it will
        % produce an error
        if locs
            allpeaks.sharpness(i,j)  = peaks(1);
            allpeaks.prominence(i,j) = prominence(1);
            % if the user specified index2time_constant we can provide
            % information regarding peak in time
            if (isfield(us_spec, 'index2time_constant'))
                allpeaks.times(i,j) = windowrange(i, 1) + (locs(1) * us_spec.index2time_constant);
                allpeaks.width(i,j) = width(1) * us_spec.index2time_constant;
            % if the user specified index2distance_constant we can provide
            % information regarding peak in distance
            elseif (isfield(us_spec, 'index2distance_constant'))
                allpeaks.locations(i,j) = windowrange(i, 1) + (locs(1) * us_spec.index2distance_constant);
                allpeaks.width(i,j)     = width(1)  * us_spec.index2distance_constant;
            end
        end

        %fprintf('frame: %d probe: %d\n', j, i);

    % end loop transducer
    end

    t_dsp(j) = toc;

% end loop frame
end

% put indicator to terminal
fprintf("DSP is finished, average Peak Detection/timeframe %.4f seconds, overall time %.4f seconds\n", ...
        mean(t_dsp), sum(t_dsp));
% close the progress bar
close(progress_bar);

% turn on the warning again
warning('on', 'signal:findpeaks:largeMinPeakHeight');

% end function
end

