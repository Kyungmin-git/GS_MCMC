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
