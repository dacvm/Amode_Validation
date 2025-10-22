function f = myfunc_fml_2a(alpha, timeGT, layer_fat, layer_muscle, layer_tenlig)
    % depthGT (mm), timeGT (mus)
    vsound_fat            = 1440.2; %m/s
    vsound_fat_reduced    = vsound_fat - (alpha*vsound_fat);
    vsound_muscle         = 1588.4; %m/s
    vsound_muscle_reduced = vsound_muscle - (alpha*vsound_muscle);
    vsound_tenlig         = 1750.0; %m/s
    vsound_tenlig_reduced = vsound_tenlig - (alpha*vsound_tenlig);


    time_fat    = ((layer_fat*1e-3)/vsound_fat_reduced) * 1e6;
    time_muscle = ((layer_muscle*1e-3)/vsound_muscle_reduced) * 1e6;
    time_tenlig = ((layer_tenlig*1e-3)/vsound_tenlig_reduced) * 1e6;

    time_est  = 2*(time_fat+time_muscle+time_tenlig);
    f = mean(abs(time_est-timeGT));
end