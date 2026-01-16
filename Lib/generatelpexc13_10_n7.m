function [candidatet0, candidateqn_n, gastype, logJacobian, logPropRatio] = generatelpexc13_10_n7(currentt0, currentqn_n)
%GENERATELPEXC13_10_N7
% Propose excitation parameters (t0, qn_n) using 4 move types:
% 1) jitter event times
% 2) jitter weights (then renormalize)
% 3) birth: add 1 impulse with Beta augmentation + time uniform
% 4) death: remove 1 impulse uniformly + renormalize
%
% Returns:
%   logJacobian   : log |J| for RJ mapping (0 for fixed-dimension moves)
%   logPropRatio  : log(q_reverse/q_forward) for MH acceptance

global tau gastype

Nt = numel(currentt0);

% --- Choose move type with guardrails on Nt
gastype = randsample(4,1);
if Nt < 2
    gastype = 3;                 % must birth
elseif Nt < 101
    gastype = randsample([1,2,3],1);
elseif Nt > 1999
    gastype = randsample([1,2,4],1);
end

% --- Tunable proposal scales
Stepsize_frac = 0.02*abs(randn());     % fraction of impulses to perturb
sigma_t_jit   = 0.04;                  % s, time jitter std
sigma_w_jit   = 0.00015;               % weight jitter std
alpha         = 1;                     % Beta/Dirichlet concentration

% For proposal ratio, you need p_B and p_D (probabilities of choosing birth/death)
% Here: reflect your selection logic approximately:
p_B = 1/4; 
p_D = 1/4;


logJacobian  = 0;
logPropRatio = 0;

switch gastype

    case 1  % ----- time jitter -----
        Nd = max(1, round(Stepsize_frac * Nt));
        ix = randperm(Nt, Nd);

        cand_t0 = currentt0;
        cand_t0(ix) = cand_t0(ix) + randn(1,Nd)*sigma_t_jit;

        % reflect into [0, tau]
        cand_t0 = abs(cand_t0);
        over = cand_t0 > tau;
        cand_t0(over) = 2*tau - cand_t0(over);

        candidatet0   = cand_t0;
        candidateqn_n = currentqn_n;

        % sort pairs
        [candidatet0, candidateqn_n] = sort_pair(candidatet0, candidateqn_n);

    case 2  % ----- weight jitter -----
        Nd = max(1, round(Stepsize_frac * Nt));
        ix = randperm(Nt, Nd);

        cand_q = currentqn_n;
        cand_q(ix) = cand_q(ix) + randn(1,Nd)*sigma_w_jit;
        cand_q = abs(cand_q);

        s = sum(cand_q);
        if s == 0
            cand_q = ones(1,Nt)/Nt;
        else
            cand_q = cand_q/s;
        end

        candidatet0   = currentt0;
        candidateqn_n = cand_q;

        % This is not RJ, so logJacobian/logPropRatio remain 0.
        [candidatet0, candidateqn_n] = sort_pair(candidatet0, candidateqn_n);

    case 3  % ----- birth (add exactly 1 impulse) -----
        % Forward proposal:
        %   t_new ~ Uniform(0,tau)
        %   x ~ Beta(alpha, Nt*alpha)
        birthn = 1;

        t_new = rand(1,birthn) * tau;

        % x in (0,1)
        x = betarnd(alpha, Nt*alpha);

        % scale old weights and append x
        cand_q_old = currentqn_n * (1 - x);
        cand_q_new = x;

        candidatet0   = [currentt0, t_new];
        candidateqn_n = [cand_q_old, cand_q_new];

        % sort pairs
        [candidatet0, candidateqn_n] = sort_pair(candidatet0, candidateqn_n);

        % Jacobian for mapping w'=(1-x)w : |J| = (1-x)^Nt
        logJacobian = Nt * log(max(1 - x, realmin));

        % Proposal ratio log(q_rev/q_fwd)
        % q_fwd = p_B * (1/tau) * BetaPDF(x; alpha, Nt*alpha)
        % q_rev = p_D * (1/(Nt+1))    [delete uniformly]
        log_q_fwd = log(p_B) + log(1/tau) + log_betapdf(x, alpha, Nt*alpha);
        log_q_rev = log(p_D) + log(1/(Nt+1));
        logPropRatio = log_q_rev - log_q_fwd;

    case 4  % ----- death (remove exactly 1 impulse) -----
        deathn = 1;

        % choose one to delete uniformly
        del_idx = randsample(Nt, deathn);
        w_del = currentqn_n(del_idx);

        % remove and renormalize
        cand_t0 = currentt0;
        cand_q  = currentqn_n;

        cand_t0(del_idx) = [];
        cand_q(del_idx)  = [];

        s = sum(cand_q);
        if s == 0
            cand_q = ones(1,numel(cand_q))/numel(cand_q);
        else
            cand_q = cand_q/s;
        end

        candidatet0   = cand_t0;
        candidateqn_n = cand_q;

        [candidatet0, candidateqn_n] = sort_pair(candidatet0, candidateqn_n);

        % Jacobian for renormalization by (1-w_del): |J| = (1/(1-w_del))^(Nt-1)
        logJacobian = (Nt-1) * log(1/max(1 - w_del, realmin));

        % Proposal ratio log(q_rev/q_fwd)
        % Forward death: q_fwd = p_D * (1/Nt)
        % Reverse birth from Nt-1 state:
        %   t_new uniform => 1/tau
        %   x ~ Beta(alpha, (Nt-1)*alpha) evaluated at x = w_del
        log_q_fwd = log(p_D) + log(1/Nt);
        log_q_rev = log(p_B) + log(1/tau) + log_betapdf(w_del, alpha, (Nt-1)*alpha);
        logPropRatio = log_q_rev - log_q_fwd;

end

end

% ----------------------- helpers -----------------------------------------
function [t_sorted, q_sorted] = sort_pair(t, q)
    M = [t(:), q(:)];
    M = sortrows(M, 1);
    t_sorted = M(:,1)'; 
    q_sorted = M(:,2)';
end

function lp = log_betapdf(x, a, b)
    % stable log Beta pdf on (0,1)
    x = min(max(x, realmin), 1-realmin);
    lp = (a-1)*log(x) + (b-1)*log(1-x) - betaln(a,b);
end
