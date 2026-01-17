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
