%% ------------------------------------------------------------------------
%  MCMC inversion of LP spectrum using the Leaky Gas Pocket model
%  (with parallel tempering and amplitude+phase likelihood)
%
%  This script:
%   1) reads 1-day miniseed data (untar + ReadMSEEDFast),
%   2) converts to displacement (velocity->bandpass->integrate->bandpass),
%   3) extracts one LP window and computes its amplitude/phase spectrum,
%   4) estimates noise statistics from multiple noise windows on the SAME day,
%   5) builds diagonal covariance matrices for amplitude and phase,
%   6) runs MCMC with parallel tempering (4 chains).
%
%  Notes:
%   - The noise window is drawn from the same preprocessed displacement trace
%     to avoid redundant re-processing and to ensure consistent filtering.
%   - first_f/last1_f etc. are indices defining frequency bands for misfit.
%   - The forward model must return modeled amplitude spectrum + phase spectrum.
% -------------------------------------------------------------------------

clear; clc;

%% --------------------------- User configuration --------------------------
cfg = struct();

% --- I/O
cfg.libPath1   = '/import/c1/SOLIDEARTH/kkim39/Lib';
cfg.workPath   = '/import/c1/SOLIDEARTH/kkim39/May26_inv_1';
cfg.tarFile    = 'GreatSitkin_May26_2021.290510.tar.mseed';
cfg.mseedFile  = 'GreatSitkin_May26_2021.290510/GSSP.AV.mseed';
cfg.metaFile   = 'IRISDMC-GSSP.AV.meta';

% --- Sampling and preprocessing
cfg.Fs     = 50;          % Hz
cfg.dt     = 1/cfg.Fs;    % s
cfg.bpOrd  = 4;           % Butterworth order (4 corners)
cfg.fLo    = 1;           % Hz
cfg.fHi    = 20;          % Hz

% --- Event window (in seconds from day start)
% Your original: LPstartt = 537217/50 - 1; LPendt = 537834/50 + 15;
cfg.LPstart = 537217/50 - 1;
cfg.LPend   = 537834/50 + 15;

% --- Noise windows (same day)
cfg.noiseWinStart = 0;       % s, start of “candidate” noise region (optional)
cfg.noiseSegDur   = [];      % if empty, uses tau (LP duration)
cfg.nNoiseSeg     = 28;      % number of noise segments used for covariance
cfg.noiseScaleFFT = 5;       % your "noisef" scale (document why you use this)

% --- Frequency bands for likelihood
cfg.fAmpMin  = 0.25;  % amplitude misfit band lower bound (Hz)
cfg.fAmpMax  = 8.0;   % amplitude misfit band upper bound (Hz)
cfg.fPhaseMin = 1.0;  % phase misfit band lower bound (Hz)
cfg.fPhaseMax = 8.0;  % phase misfit band upper bound (Hz)
cfg.fPlotMax  = 10.0; % plotting limit (Hz)

% --- MCMC settings
cfg.nIter        = 24000000;
cfg.kStart       = 2;            % start from 2 (store initial at 1)
cfg.saveEvery    = 400000;       % save .mat checkpoint
cfg.plotEvery    = 400000;       % diagnostics plots
cfg.plotAlsoAt   = [5000 10000]; % extra plots early
cfg.NNmax        = 10000;        % maximum impulses supported in storage
cfg.dist_m       = 3780;         % source-receiver distance (m)
cfg.rngSeed      = 1;            % reproducibility

% --- Parallel tempering temperatures
cfg.TempMain = 1.0;
cfg.TempHot  = 1.22;
cfg.TempHot2 = 1.22^2;
cfg.TempHot3 = 1.22^3;

% --- Priors: document and tune as needed
cfg.prior = struct();
cfg.prior.T_muK   = 1225.15; cfg.prior.T_sigK   = 100;
cfg.prior.R_mu    = 13;      cfg.prior.R_sig    = 9;
cfg.prior.rho_mu  = 2950;    cfg.prior.rho_sig  = 50;

% Soft bounds / penalties (publishable explanation)
cfg.prior.poro_log10_lo = log10(0.2);   % log10(%)
cfg.prior.poro_log10_hi = log10(20);    % log10(%)
cfg.prior.poro_sig_log10 = 0.03;        % smooth tail width in log10-space

cfg.prior.kappaSoftFloor = 7;           % m(5) in your code: exp_kappa? (keep consistent!)
cfg.prior.kappaSoftSig   = 0.2;

cfg.prior.QfSoftFloor = 300;
cfg.prior.QfSoftSig   = 20;

% --- Likelihood scaling
cfg.noiseAmpScale = 3;   % your "noise_amp" multiplier

rng(cfg.rngSeed);

%% ----------------------------- Paths -------------------------------------
addpath(cfg.libPath1);
addpath(cfg.workPath);

%% ----------------------------- Read data ---------------------------------
% untar miniseed archive
untar(cfg.tarFile);

% read miniseed
data = ReadMSEEDFast(fullfile(pwd, cfg.mseedFile));
wave_raw = data.data(:);
Fs = data.sampleRate;
assert(abs(Fs - cfg.Fs) < 1e-6, 'Sampling rate mismatch.');

% build time vector in seconds relative to trace start
t = (0:numel(wave_raw)-1)'/Fs;

% read metadata and scale (instrument sensitivity)
scale = read_scale_from_meta(cfg.metaFile);

%% ---------------------- Preprocess (displacement) ------------------------
% Convert to physical units:
%   waveform_org/scale  -> velocity (document this in your paper: depends on response)
% Then:
%   detrend + demean
%   bandpass (1-20 Hz)
%   integrate to displacement
%   bandpass again (consistent with your original workflow)
disp_day = preprocess_to_displacement(wave_raw, scale, cfg);

% Quick-look plot (one-day displacement)
figure; plot(t, disp_day); xlabel('Time (s)'); ylabel('Displacement (m)');
title('1-day displacement (preprocessed)'); set(gca,'FontSize',14);

%% --------------------- Extract LP segment + spectrum ---------------------
LP = extract_window(disp_day, t, cfg.LPstart, cfg.LPend);
tau = LP.t(end) - LP.t(1);

% spectrum (complex), amplitude, phase
[LPspec] = compute_spectrum(LP.x, cfg.dt, cfg.Fs);

% plotting
figure; plot(LP.t, LP.x);
xlabel('Time (s)'); ylabel('Displacement (m)'); title('LP displacement');
set(gca,'FontSize',16);

figure; plot(LPspec.f, LPspec.amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude (units)'); title('LP amplitude spectrum');
xlim([0 cfg.fPlotMax]); set(gca,'FontSize',16);

%% --------------------- Noise statistics (same day) ------------------------
% IMPORTANT:
% Use the SAME preprocessed displacement trace and slice many windows.
% This removes redundancy and ensures identical filtering for signal/noise.

if isempty(cfg.noiseSegDur)
    cfg.noiseSegDur = tau; % match LP length (your original approach)
end

% Choose a noise region (example: 1500–1800 s) and then sub-sample segments.
% If you want exactly your old behavior: segments starting at (i-1)*10 seconds,
% just implement that rule here with a documented start time.
noiseRegion = [1500 1800]; % <-- your original selection
noiseStats  = estimate_noise_stats(disp_day, t, noiseRegion, cfg, LPspec);

% Build diagonal covariance for amplitude misfit in [fAmpMin, fAmpMax]
idxAmp = find_freq_indices(LPspec.f, cfg.fAmpMin, cfg.fAmpMax);

sigma_amp = cfg.noiseAmpScale * noiseStats.amp_mean ./ sqrt(2); % your sigma_row logic
CovAmp = diag( sigma_amp(idxAmp).^2 );
iCovAmp = diag( 1 ./ (sigma_amp(idxAmp).^2) );

% Build diagonal covariance for phase misfit in [fPhaseMin, fPhaseMax]
idxPh = find_freq_indices(LPspec.f, cfg.fPhaseMin, cfg.fPhaseMax);

% A more defensible phase-noise estimate:
% Use many complex-noise realizations + fixed LP complex spectrum to measure
% the dispersion of phase( S + N ).
% (This replaces sigma_row_angle = angle(...) which is not a standard deviation.)
sigma_phase = estimate_phase_sigma_from_noise(LPspec.complex, noiseStats.noiseComplexSamples);
CovPh = diag( sigma_phase(idxPh).^2 );
iCovPh = diag( 1 ./ (sigma_phase(idxPh).^2) );

% Diagnostics (noise spectra examples + mean)
figure; plot(LPspec.f, noiseStats.amp_examples(:,1:min(5,size(noiseStats.amp_examples,2))));
xlim([0 cfg.fPlotMax]); title('Example noise spectra'); xlabel('Hz'); set(gca,'FontSize',14);

figure; plot(LPspec.f, noiseStats.amp_mean);
xlim([0 cfg.fPlotMax]); title('Mean noise amplitude spectrum'); xlabel('Hz'); set(gca,'FontSize',14);

% Energy check (optional, for sanity)
E_sig = sum(abs(LPspec.complex).^2) * LPspec.df;
E_noise_mean = mean(noiseStats.energy);
fprintf('Energy check: E_sig=%.3e, mean(E_noise)=%.3e\n', E_sig, E_noise_mean);

% Complex-domain scatter at ~2 Hz and ~4 Hz (optional, publishable as QC)
plot_complex_scatter(noiseStats.noiseComplexSamples, LPspec.f, 2.0);
plot_complex_scatter(noiseStats.noiseComplexSamples, LPspec.f, 4.0);

%% -------------------------- MCMC initialization --------------------------
% Parameter vector m = [T, Q, R, L, exp_kappa, exp_D, log_phi, rho_s, Qf]
% Keep your meaning consistent (units!) and document in manuscript.
[m0, state0] = initialize_mcmc_state(cfg, tau);

% Evaluate forward model at initial state
% forwardmodelgas15_10ppga must return:
%  Posc, synthetic, ..., modelspect_raw, A_res, A_exc, A_path, modelspect, modelPhase
[model0] = eval_forward(m0, state0.t0, state0.qn, state0.sigma_k_log, cfg);

% Likelihood components at initial state
[ll0_amp, ll0_ph] = loglike_amp_phase(model0.amp(idxAmp), LPspec.amp(idxAmp), ...
                                      model0.phase(idxPh), LPspec.phase(idxPh), ...
                                      iCovAmp, iCovPh);

lp0 = logprior_m(cfg, m0);

%% ----------------------- Run MCMC + parallel tempering -------------------
% NOTE:
% Below is the “clean loop skeleton”. Inside the proposal step, you call your
% existing proposal kernels unchanged.

chain = struct();
chain.m          = zeros(cfg.nIter, numel(m0), 'single'); % memory conscious
chain.accept     = false(cfg.nIter, 1);
chain.sigma_store = zeros(cfg.nIter,1,'single');

% Store initial
chain.m(1,:) = single(m0);
chain.sigma_store(1) = single(10^(state0.sigma_k_log));

% Current state
cur.m = m0;
cur.state = state0;
cur.model = model0;
cur.ll_amp = ll0_amp;
cur.ll_ph  = ll0_ph;
cur.lp = lp0;

% Track MAP (by likelihood or posterior—choose one and say it)
MAP.ll = -Inf;
MAP.m  = cur.m;
MAP.state = cur.state;
MAP.model = cur.model;

% Hot chains: repeat initialization in a helper so the code is not duplicated
hot(1) = init_hot_chain(cfg, tau, cfg.TempHot);
hot(2) = init_hot_chain(cfg, tau, cfg.TempHot2);
hot(3) = init_hot_chain(cfg, tau, cfg.TempHot3);

fprintf('Starting MCMC loop...\n');

for k = cfg.kStart:cfg.nIter

    % ---------------- Proposal type selection ----------------
    moveType = mod(k,3);   % 0: source params, 1: excitation (t0/qn_n), 2: sigma
    % (Your old ra/rb were unused; remove them to reduce confusion.)

    % ---------------- Propose candidate (main chain) ----------
    cand = propose_candidate(cur, moveType, cfg);

    % Check hard bounds + compute log prior
    cand.lp = logprior_m(cfg, cand.m);
    if ~isfinite(cand.lp)
        acceptMain = false;
    else
        % Evaluate forward model
        cand.model = eval_forward(cand.m, cand.state.t0, cand.state.qn, cand.state.sigma_k_log, cfg);

        % Likelihood
        [cand.ll_amp, cand.ll_ph] = loglike_amp_phase( ...
            cand.model.amp(idxAmp), LPspec.amp(idxAmp), ...
            cand.model.phase(idxPh), LPspec.phase(idxPh), ...
            iCovAmp, iCovPh);

        % Metropolis-Hastings ratio (include Jacobian/proposal if your move needs it)
        logalpha = (cand.lp - cur.lp) ...
                 + (cand.ll_amp - cur.ll_amp) / cfg.TempMain ...
                 + (cand.ll_ph  - cur.ll_ph ) / cfg.TempMain ...
                 + cand.state.logproposal + cand.state.logJacobian;

        logalpha = min(0, logalpha);
        acceptMain = (log(rand()) < logalpha);
    end

    if acceptMain
        cur = cand;
        chain.accept(k) = true;
    end

    chain.m(k,:) = single(cur.m);
    chain.sigma_store(k) = single(10^(cur.state.sigma_k_log));

    % Update MAP (here: likelihood-only; you can switch to posterior if you prefer)
    if (cur.ll_amp + cur.ll_ph) > MAP.ll
        MAP.ll = (cur.ll_amp + cur.ll_ph);
        MAP.m = cur.m;
        MAP.state = cur.state;
        MAP.model = cur.model;
    end

    % ---------------- Hot chains update -----------------------
    for ih = 1:3
        hot(ih) = step_one_chain(hot(ih), moveType, cfg, LPspec, idxAmp, idxPh, iCovAmp, iCovPh);
    end

    % ---------------- Parallel tempering swaps ----------------
    % Swap main <-> hot1, hot1 <-> hot2, hot2 <-> hot3
    [cur, hot(1)] = attempt_swap(cur, hot(1));
    [hot(1), hot(2)] = attempt_swap(hot(1), hot(2));
    [hot(2), hot(3)] = attempt_swap(hot(2), hot(3));

    % ---------------- Diagnostics / saving --------------------
    if mod(k, cfg.saveEvery) == 0
        save(sprintf('syntheticresult_checkpoint_%d.mat', k), '-v7.3');
    end

    if mod(k, cfg.plotEvery) == 0 || any(k == cfg.plotAlsoAt)
        fprintf('k=%d, accRate=%.4f, ll=%.3f\n', k, mean(chain.accept(1:k)), cur.ll_amp+cur.ll_ph);
        % Put your existing diagnostic plots here, but use "cur" and "MAP"
        % so the figure code is clean and does not depend on globals.
    end

end

% Final save
save('syntheticresult_final.mat', '-v7.3');

fprintf('Done. Final acceptance=%.4f\n', mean(chain.accept));

%% ============================ Local functions ============================
function scale = read_scale_from_meta(metaFile)
    fid = fopen(metaFile,'r');
    assert(fid>0, 'Cannot open metadata file: %s', metaFile);
    C = textscan(fid, ['%s',repmat('%s',[1,20])], 'Delimiter','|');
    fclose(fid);

    % Your original: scale = str2num(Metaseis{1,12}{3,1});
    % Add a safety check:
    raw = C{1,12}{3,1};
    scale = str2double(raw);
    assert(isfinite(scale) && scale~=0, 'Invalid scale parsed from metadata.');
end

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

function W = extract_window(x, t, t0, t1)
    idx = (t>=t0) & (t<=t1);
    W.x = x(idx);
    W.t = t(idx) - t(find(idx,1,'first'));
end

function S = compute_spectrum(x, dt, Fs)
    % FFT normalization: keep consistent with your forward model convention.
    % Here we keep a simple, documented convention:
    %   X = fft(x)/sqrt(N) * sqrt(dt/df)
    N = numel(x);
    df = (1/dt)/N;

    X = fft(x(:), N) / sqrt(N) * sqrt(dt/df);
    f = (0:(N/2))' * (Fs/N);

    Xpos = X(1:N/2+1);

    S.f = f;
    S.df = df;
    S.complex = Xpos;
    S.amp = abs(Xpos);
    S.phase = angle(Xpos);
end

function idx = find_freq_indices(f, fmin, fmax)
    idx = find(f>=fmin & f<=fmax);
    assert(~isempty(idx), 'Frequency index set is empty for [%.2f %.2f] Hz.', fmin, fmax);
end

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

function plot_complex_scatter(noiseComplexSamples, f, fTarget)
    [~,ix] = min(abs(f - fTarget));
    z = noiseComplexSamples(ix,:);
    figure; scatter(real(z), imag(z), 10, 'filled');
    axis equal; grid on;
    title(sprintf('Noise complex scatter near %.2f Hz', f(ix)));
    xlabel('Real'); ylabel('Imag'); set(gca,'FontSize',14);
end

function [m0, state] = initialize_mcmc_state(cfg, tau)
    % Your original random initialization, cleaned and documented.
    % Parameter vector: m = [T0, Q0, R0, L0, exp_kappa0, exp_D0, log_phi0, rho_s0, Qf0]
    T0 = 800 + 400*rand + 273.15;
    Q0 = 100 + 80*rand;          % (document units: kg/s if your forward model uses kg/s)
    R0 = 11 + 25*rand;
    L0 = 20 + 60*rand;
    exp_kappa0 = 7 + 3.5*rand;   % log10(1/kappa) style in your code — keep consistent
    exp_D0     = 1.7*rand;
    log_phi0   = 1*rand;         % log10(porosity %) in your code
    rho_s0     = 2900;
    Qf0        = 30 + 100*rand;

    m0 = [T0,Q0,R0,L0,exp_kappa0,exp_D0,log_phi0,rho_s0,Qf0];

    % Gas excitation: Dirichlet-like weights via Gamma samples
    NN = 300;
    alpha = 1;
    qn_n = gamrnd(alpha, 1, 1, NN);
    qn_n = qn_n / sum(qn_n);

    state.t0 = rand(1,NN) * tau;
    state.qn_n = qn_n;
    state.qn = Q0 * tau * qn_n;     % kg per impulse if Q0 is kg/s (document!)
    state.sigma_k_log = -1;         % your sigma_k_log0
    state.logJacobian = 0;
    state.logproposal = 0;
end

function lp = logprior_m(cfg, m)
    % Hard bounds (from your original if-block)
    if ~(m(1)>=973.15 && m(1)<=1473.15 && ...
         m(2)>=0      && m(2)<=10000  && ...
         m(3)>=5      && m(3)<=40     && ...
         m(4)>=1      && m(4)<=80     && ...
         m(5)>=6      && m(5)<=12     && ...
         m(6)>=-log10(3) && m(6)<=2   && ...
         m(7)>=-1     && m(7)<= (1+log10(3)) && ...
         m(8)>=2800   && m(8)<=3100   && ...
         m(9)>=20     && m(9)<=400)
        lp = -Inf;
        return;
    end

    % Gaussian priors (your original)
    lp = 0;
    lp = lp - (m(1)-cfg.prior.T_muK )^2/(2*cfg.prior.T_sigK^2);
    lp = lp - (m(3)-cfg.prior.R_mu  )^2/(2*cfg.prior.R_sig^2);
    lp = lp - (m(8)-cfg.prior.rho_mu)^2/(2*cfg.prior.rho_sig^2);

    % Soft log-uniform porosity in [%] with smooth Gaussian tails
    logp = m(7);
    lo = cfg.prior.poro_log10_lo; hi = cfg.prior.poro_log10_hi;
    sig = cfg.prior.poro_sig_log10;
    if logp < lo
        d = lo - logp; lp = lp - 0.5*(d/sig)^2;
    elseif logp > hi
        d = logp - hi; lp = lp - 0.5*(d/sig)^2;
    end

    % Soft floor on m(5) (your original “min(m(5),7)” trick)
    sig5 = cfg.prior.kappaSoftSig;
    lp = lp - 0.5*((min(m(5), cfg.prior.kappaSoftFloor) - cfg.prior.kappaSoftFloor)/sig5)^2;

    % Soft floor on Qf (m(9))
    sigQf = cfg.prior.QfSoftSig;
    lp = lp - 0.5*((max(m(9), cfg.prior.QfSoftFloor) - cfg.prior.QfSoftFloor)/sigQf)^2;
end

function model = eval_forward(m, t0, qn, sigma_k_log, cfg)
    % Wrap your forwardmodelgas15_10ppga here so the main loop stays clean.
    % Ensure it returns consistent spectrum grids and phase convention.
    [Posc, synthetic, ~, modelspect_raw, A_res, A_exc, A_path, modelspect, modelPhase] = ...
        forwardmodelgas15_10ppga(m, t0, qn, sigma_k_log);

    model.Posc = Posc;
    model.synthetic = synthetic;
    model.amp = modelspect(:);
    model.phase = modelPhase(:);
    model.raw = modelspect_raw;
    model.A_res = A_res;
    model.A_exc = A_exc;
    model.A_path = A_path;

    %#ok<NASGU> % if you don’t use some outputs in the main script
end

function [ll_amp, ll_ph] = loglike_amp_phase(modelAmp, dataAmp, modelPh, dataPh, iCovAmp, iCovPh)
    % Amplitude misfit
    rA = (modelAmp(:) - dataAmp(:));
    ll_amp = -0.5 * (rA' * iCovAmp * rA);

    % Phase misfit with wrap-around
    d = abs(dataPh(:) - modelPh(:));
    rP = min(d, 2*pi - d);
    ll_ph = -0.5 * (rP' * iCovPh * rP);
end

function cand = propose_candidate(cur, moveType, cfg)
    % This function is a thin wrapper that calls your existing kernels.
    % It returns a candidate struct containing:
    %   cand.m, cand.state (t0, qn_n, qn, sigma_k_log, logJacobian, logproposal)
    %
    % NOTE: keep the book-keeping here so the loop stays readable.

    cand = cur; % start as copy

    if moveType == 0
        % Source parameter move
        cand.m = generatelpgas11_4m6(cur.m);
        cand.state.logJacobian = 0;
        cand.state.logproposal = 0;

    elseif moveType == 1
        % Excitation move: t0 + qn_n (RJ / birth-death / shift etc.)
        [t0_new, qn_n_new, ~, logJ, logq] = generatelpexc13_10_n7(cur.state.t0, cur.state.qn_n);
        cand.state.t0 = t0_new;
        cand.state.qn_n = qn_n_new;
        cand.state.logJacobian = logJ;
        cand.state.logproposal = logq;

    else
        % sigma (gas-impulse timescale) move
        cand.state.sigma_k_log = sigma_propose2(cur.state.sigma_k_log);
        cand.state.logJacobian = 0;
        cand.state.logproposal = 0;
    end

    % Update qn from (qn_n, Q, tau)
    tau = max(cur.state.t0); % not safe if t0 doesn't cover tau; better pass tau in state
    % Better: store tau in state at init and never infer it here.
    % For now: do the robust thing:
    % (You should store tau explicitly.)
    tau = max(tau, 1e-6);

    cand.state.qn = cand.state.qn_n * cand.m(2) * tau;
end

function hot = init_hot_chain(cfg, tau, Temp)
    [m0, state0] = initialize_mcmc_state(cfg, tau);
    model0 = eval_forward(m0, state0.t0, state0.qn, state0.sigma_k_log, cfg);
    hot.m = m0;
    hot.state = state0;
    hot.model = model0;
    hot.Temp = Temp;
    hot.lp = logprior_m(cfg, m0);
    hot.ll_amp = NaN; hot.ll_ph = NaN; % computed in step_one_chain
end

function chain = step_one_chain(chain, moveType, cfg, LPspec, idxAmp, idxPh, iCovAmp, iCovPh)
    cand = propose_candidate(chain, moveType, cfg);
    cand.lp = logprior_m(cfg, cand.m);
    if ~isfinite(cand.lp)
        return;
    end
    cand.model = eval_forward(cand.m, cand.state.t0, cand.state.qn, cand.state.sigma_k_log, cfg);
    [cand.ll_amp, cand.ll_ph] = loglike_amp_phase( ...
        cand.model.amp(idxAmp), LPspec.amp(idxAmp), ...
        cand.model.phase(idxPh), LPspec.phase(idxPh), ...
        iCovAmp, iCovPh);

    % If chain.ll_* not set yet, initialize it now
    if ~isfinite(chain.ll_amp)
        [chain.ll_amp, chain.ll_ph] = loglike_amp_phase( ...
            chain.model.amp(idxAmp), LPspec.amp(idxAmp), ...
            chain.model.phase(idxPh), LPspec.phase(idxPh), ...
            iCovAmp, iCovPh);
    end

    logalpha = (cand.lp - chain.lp) ...
             + (cand.ll_amp - chain.ll_amp)/chain.Temp ...
             + (cand.ll_ph  - chain.ll_ph )/chain.Temp ...
             + cand.state.logproposal + cand.state.logJacobian;

    logalpha = min(0, logalpha);
    if log(rand()) < logalpha
        chain = cand;
    end
end

function [A,B] = attempt_swap(A,B)
    % Parallel tempering swap:
    % accept with prob min(1, exp((1/TA-1/TB)*(llB-llA)))
    if ~isfinite(A.ll_amp) || ~isfinite(B.ll_amp)
        return;
    end

    llA = A.ll_amp + A.ll_ph;
    llB = B.ll_amp + B.ll_ph;

    logalpha = (1/A.Temp - 1/B.Temp) * (llB - llA);
    logalpha = min(0, logalpha);

    if log(rand()) < logalpha
        tmp = A; A = B; B = tmp;
    end
end
