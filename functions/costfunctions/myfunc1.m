function f = myfunc1(alpha, timeGT, depthGT)
    % depthGT (mm), timeGT (mus)
    vsound_softtissue = 1540; %m/s
    vsound_reducted   = vsound_softtissue - (alpha*vsound_softtissue);

    time_est = 2 * ( (depthGT*1e-3)/vsound_reducted) * 1e6;
    f = mean(abs(time_est-timeGT));
end