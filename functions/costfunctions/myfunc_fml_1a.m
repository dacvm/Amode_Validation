function f = myfunc_fml_1a(alpha, timeGT, layer_fat, layer_muscle)
    % depthGT (mm), timeGT (mus)
    vsound_fat            = 1440.2; %m/s
    vsound_fat_reduced    = vsound_fat - (alpha*vsound_fat);
    vsound_muscle         = 1588.4; %m/s
    vsound_muscle_reduced = vsound_muscle - (alpha*vsound_muscle);

    time_fat    = ((layer_fat*1e-3)/vsound_fat_reduced) * 1e6;
    time_muscle = ((layer_muscle*1e-3)/vsound_muscle_reduced) * 1e6;

    time_est  = 2*(time_fat+time_muscle);
    f = mean(abs(time_est-timeGT));
end