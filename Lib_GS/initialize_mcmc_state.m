function [m0, state] = initialize_mcmc_state(cfg, tau)
    % Your original random initialization, cleaned and documented.
    % Parameter vector: m = [T0, Q0, R0, L0, exp_kappa0, exp_D0, log_phi0, rho_s0, Qf0]
    T0 = 800 + 400*rand + 273.15;
    Q0 = 200 + 180*rand;          % (document units: kg/s if your forward model uses kg/s)
    R0 = 11 + 25*rand;
    L0 = 40 + 30*rand;
    exp_kappa0 = 7 + 3.5*rand;   % log10(1/kappa) style in your code — keep consistent
    exp_D0     = 1.7*rand;
    log_phi0   = 1*rand;         % log10(porosity %) in your code
    rho_s0     = 2900;
    Qf0        = 150 + 100*rand;

    m0 = [T0,Q0,R0,L0,exp_kappa0,exp_D0,log_phi0,rho_s0,Qf0];

    % Gas excitation: Dirichlet-like weights via Gamma samples
    NN = 300;
    alpha = 1;
    qn_n = gamrnd(alpha, 1, 1, NN);
    qn_n = qn_n / sum(qn_n);

    state.t0 = rand(1,NN) * tau;
    state.qn_n = qn_n;
    state.qn = Q0 * tau * qn_n;     % kg per impulse if Q0 is kg/s (document!)
    state.sigma_k_log = -1.5;         % your sigma_k_log0
    state.logJacobian = 0;
    state.logproposal = 0;
    state.tau = tau;   % store tau explicitly

end
