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

 clc;

%% --------------------------- User configuration --------------------------
cfg = struct();

% --- I/O
cfg.libPath1   = '/Users/kyuninkim/Downloads/GS_MCMC_code1/Lib_GS';
cfg.workPath   = '/Users/kyuninkim/Downloads/GS_MCMC_code1/MCMC_cluster18';
cfg.datapath   = '/Users/kyuninkim/Downloads/GS_MCMC_code1/Data_file';
cfg.tarFile    = 'GreatSitkin_August08_2021_0000.tar.mseed';
cfg.mseedFile  = 'GreatSitkin_August08_0000.240701/GSSP.AV.mseed';
cfg.metaFile   = 'IRISDMC-GSSP.AV.meta';

% --- Sampling and preprocessing
cfg.Fs     = 50;          % Hz
cfg.dt     = 1/cfg.Fs;    % s
cfg.bpOrd  = 4;           % Butterworth order (4 corners)
cfg.fLo    = 1;           % Hz
cfg.fHi    = 20;          % Hz

% --- Event window (in seconds from day start)
% Your original: LPstartt = 537217/50 - 1; LPendt = 537834/50 + 15;
cfg.LPstart = 3864115/50 - 1;
cfg.LPend   = 3864562/50 + 15;

% --- Noise windows (same day)
cfg.noiseWinStart = 0;       % s, start of “candidate” noise region (optional)
cfg.noiseSegDur   = [];      % if empty, uses tau (LP duration)
cfg.nNoiseSeg     = 28;      % number of noise segments used for covariance


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
global distance stepsize stepsize2
distance =cfg.dist_m;
cfg.rngSeed      = 1;            % reproducibility
stepsize2=0.003;
stepsize=0.003*[1.2*400,60,1*35,1*70,0.5*3,0.4*2.5,1.2,1*200,1.1*100];
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
cfg.noiseAmpScale = 1;   
rng(cfg.rngSeed);
kk=3;

%% ----------------------------- Paths -------------------------------------
addpath(cfg.libPath1);
addpath(cfg.workPath);
addpath(cfg.datapath);

%% ----------------------------- Read data ---------------------------------
% untar miniseed archive
untar(cfg.tarFile);

% Find the extracted directory
d = dir(cfg.tarFile);
tarBase = erase(cfg.tarFile, '.tar.mseed');

dataDir = fullfile(pwd, tarBase);

assert(isfolder(dataDir), 'Extracted folder not found: %s', dataDir);


% read miniseed
data = ReadMSEEDFast(fullfile(pwd, cfg.mseedFile));
wave_raw = data.data(:);
Fs = data.sampleRate;
assert(abs(Fs - cfg.Fs) < 1e-6, 'Sampling rate mismatch.');

% build time vector in seconds relative to trace start
t = (0:numel(wave_raw)-1)'/Fs;

% read metadata and scale (instrument sensitivity)
scale = read_scale_from_meta(fullfile(dataDir, cfg.metaFile));

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
global tau
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
%sigma_phase = estimate_phase_sigma_from_noise(LPspec.complex, noiseStats.noiseComplexSamples);
noiseComplexScaled = cfg.noiseAmpScale * noiseStats.noiseComplexSamples;
sigma_phase = estimate_phase_sigma_from_noise(LPspec.complex, noiseComplexScaled);


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
% Keep your meaning consistent and document in manuscript.
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
% Current state
cur.m = m0;
cur.state = state0;
cur.model = model0;
cur.ll_amp = ll0_amp;
cur.ll_ph  = ll0_ph;
cur.lp = lp0;

cur.Temp = cfg.TempMain;   % <-- ADD THIS LINE

% Track MAP (by likelihood or posterior—choose one and say it)
MAP.ll = -Inf;
MAP.m  = cur.m;
MAP.state = cur.state;
MAP.model = cur.model;

% Hot chains: repeat initialization in a helper so the code is not duplicated
hot(1) = init_hot_chain(cfg, tau, cfg.TempHot);
hot(2) = init_hot_chain(cfg, tau, cfg.TempHot2);
hot(3) = init_hot_chain(cfg, tau, cfg.TempHot3);

    % allocate once before loop
chain.ll = nan(cfg.nIter,1,'single');     % total (amp+phase)
chain.ll_amp = nan(cfg.nIter,1,'single');
chain.ll_ph  = nan(cfg.nIter,1,'single');


fprintf('Starting MCMC loop...\n');
disp('cur.state fields:'); disp(fieldnames(cur.state));
disp('cur.state.sigma_k_log:'); disp(cur.state.sigma_k_log);

for k = cfg.kStart:cfg.nIter
    k
    % ---------------- Proposal type selection ----------------
    moveType = mod(k,3);   % 0: source params, 1: excitation (t0/qn_n), 2: sigma
   

    % ---------------- Propose candidate (main chain) ----------
    cand = propose_candidate(cur, moveType, cfg);

if ~isfield(cand.state,'sigma_k_log') || isempty(cand.state.sigma_k_log)
    error('sigma_k_log dropped RIGHT AFTER propose_candidate at k=%d, moveType=%d', k, moveType);
end


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

    % inside loop, after possible accept:
    chain.ll(k)     = single(cur.ll_amp + cur.ll_ph);
    chain.ll_amp(k) = single(cur.ll_amp);
    chain.ll_ph(k)  = single(cur.ll_ph);

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
        save('syntheticresult_checkpoint.mat'), '-v7.3');
    end

    if mod(k, cfg.plotEvery) == 0 || any(k == cfg.plotAlsoAt)
        fprintf('k=%d, accRate=%.4f, ll=%.3f\n', k, mean(chain.accept(1:k)), cur.ll_amp+cur.ll_ph);
        
        do_diagnostics(k, cfg, chain, cur, MAP, LP, LPspec, idxAmp, idxPh,kk);
    end

end

% Final save
save('syntheticresult_final.mat', '-v7.3');

fprintf('Done. Final acceptance=%.4f\n', mean(chain.accept));