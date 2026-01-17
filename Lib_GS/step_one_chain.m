
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
