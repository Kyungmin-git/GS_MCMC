
function cand = propose_candidate(cur, moveType, cfg)
    % Copy everything (including sigma_k_log) first
    cand = cur;

    % Always reset RJ bookkeeping unless the move sets them
    cand.state.logJacobian = 0;
    cand.state.logproposal = 0;

    if moveType == 0
        % Source parameter move (m only)
        cand.m = generatelpgas11_4m6(cur.m);

    elseif moveType == 1
        % Excitation move (t0, qn_n only)
        [t0_new, qn_n_new, ~, logJ, logq] = generatelpexc13_10_n7(cur.state.t0, cur.state.qn_n);

        cand.state.t0   = t0_new;
        cand.state.qn_n = qn_n_new;
        cand.state.logJacobian = logJ;
        cand.state.logproposal = logq;

        % IMPORTANT: do NOT touch cand.state.sigma_k_log here
        % It should carry over from cur by virtue of cand=cur above.

    else
        % sigma move (sigma only)
        cand.state.sigma_k_log = sigma_propose2(cur.state.sigma_k_log);
    end

    % ---- Required-field enforcement (debug + safety) ----
    cand = ensure_state_fields(cand);

    % Update qn using stored tau (do NOT infer tau from t0)
    tau = cand.state.tau;
    cand.state.qn = cand.state.qn_n * cand.m(2) * tau;
end