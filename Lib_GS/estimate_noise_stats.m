function stats = estimate_noise_stats(disp_day, t, noiseRegion, cfg, LPspec)
    % Extract a region, then slice multiple segments of duration cfg.noiseSegDur.
    t0 = noiseRegion(1); t1 = noiseRegion(2);
    assert(t1>t0, 'Invalid noiseRegion');

    segDur = cfg.noiseSegDur;
    Nseg = cfg.nNoiseSeg;

    % Evenly place start times within noise region (or implement your old (i-1)*10 rule)
    latestStart = t1 - segDur;
    starts = linspace(t0, latestStart, Nseg);

    % Store amplitude spectra examples and complex spectra for phase stats
    ampMat = zeros(numel(LPspec.f), Nseg);
    complexMat = zeros(numel(LPspec.f), Nseg);
    energy = zeros(1,Nseg);

    for i = 1:Nseg
        W = extract_window(disp_day, t, starts(i), starts(i)+segDur);
        Si = compute_spectrum(W.x, cfg.dt, cfg.Fs);

        % Interpolate to the LP frequency grid if lengths differ (shouldn’t if segDur matches tau)
        if numel(Si.f) ~= numel(LPspec.f) || any(abs(Si.f - LPspec.f) > 1e-12)
            ampMat(:,i) = interp1(Si.f, Si.amp, LPspec.f, 'linear', 'extrap');
            complexMat(:,i) = interp1(Si.f, Si.complex, LPspec.f, 'linear', 'extrap');
        else
            ampMat(:,i) = Si.amp;
            complexMat(:,i) = Si.complex;
        end

        energy(i) = sum(abs(complexMat(:,i)).^2) * LPspec.df;
    end

    stats.amp_examples = ampMat;
    stats.amp_mean = mean(ampMat, 2);
    stats.noiseComplexSamples = complexMat;
    stats.energy = energy;
end
