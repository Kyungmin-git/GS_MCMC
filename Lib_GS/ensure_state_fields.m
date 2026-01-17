function S = ensure_state_fields(S)
    % Ensure required state fields exist and are scalar where needed.
    if ~isfield(S,'state') || isempty(S.state)
        error('State struct missing entirely.');
    end

    if ~isfield(S.state,'sigma_k_log') || isempty(S.state.sigma_k_log)
        % If it ever happens, this tells you the bug is earlier than this call.
        error('sigma_k_log missing inside ensure_state_fields (state got overwritten somewhere).');
    end
    if ~isscalar(S.state.sigma_k_log)
        error('sigma_k_log must be scalar.');
    end

    if ~isfield(S.state,'tau') || isempty(S.state.tau) || ~isscalar(S.state.tau)
        error('tau missing from state (store tau at init and carry it).');
    end

    if ~isfield(S.state,'logJacobian') || isempty(S.state.logJacobian)
        S.state.logJacobian = 0;
    end
    if ~isfield(S.state,'logproposal') || isempty(S.state.logproposal)
        S.state.logproposal = 0;
    end
end


