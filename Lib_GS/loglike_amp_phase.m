function [ll_amp, ll_ph] = loglike_amp_phase(modelAmp, dataAmp, modelPh, dataPh, iCovAmp, iCovPh)
    % Amplitude misfit
    rA = (modelAmp(:) - dataAmp(:));
    ll_amp = -0.5 * (rA' * iCovAmp * rA);

    % Phase misfit with wrap-around
    d = abs(dataPh(:) - modelPh(:));
    rP = min(d, 2*pi - d);
    ll_ph = -0.5 * (rP' * iCovPh * rP);
end

