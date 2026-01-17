function sigma_phase = estimate_phase_sigma_from_noise(S_complex, noiseComplexSamples)
    % Estimate phase uncertainty at each frequency by Monte Carlo:
    %   phi_i = angle(S + N_j)
    % sigma_phase = circular standard deviation of phi_i around angle(S)
    %
    % This is more defensible than treating an angle() output as a “sigma”.

    phi0 = angle(S_complex);
    phiSamp = angle(S_complex + noiseComplexSamples);

    % Circular difference to phi0 in [-pi, pi]
    dphi = angle(exp(1i*(phiSamp - phi0)));

    % Circular std approximation: std of wrapped differences
    sigma_phase = std(dphi, 0, 2);

    % Avoid zeros (can blow up iCov)
    sigma_phase(sigma_phase < 1e-3) = 1e-3;
end
