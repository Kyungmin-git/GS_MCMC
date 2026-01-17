function x_disp = preprocess_to_displacement(x_raw, scale, cfg)
    % Convert to velocity-like physical units (depending on instrument response)
    x = x_raw(:) / scale;

    % Detrend + demean
    x = detrend(x);
    x = x - mean(x);

    % Bandpass (Butterworth)
    [b,a] = butter(cfg.bpOrd/2, [cfg.fLo cfg.fHi]/(cfg.Fs/2), 'bandpass');
    x = filtfilt(b,a,x);

    % Integrate to displacement
    x = cumtrapz(cfg.dt, x);

    % Bandpass again (your original workflow)
    x_disp = filtfilt(b,a,x);
end