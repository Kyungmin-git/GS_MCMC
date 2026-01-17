function do_diagnostics(kk, cfg, chain, cur, MAP, LP, LPspec, idxAmp, idxPh)
% Compatible with the updated struct-based code.
% Saves:
%   llcurrent_k.jpg
%   MCMC_source_k.jpg
%   fit_synthetic_k.jpg
%   fit_timedomain_k.jpg

    % ---------- settings ----------
    burn = 1001;
    if k < burn
        burn = 1;
    end
    mout = double(chain.m(burn:k,:));  % [T Q R L exp_kappa exp_D log_phi rho Qf]
    f = LPspec.f(:);
    tau = cur.state.tau;

    % ---------- 1) log-likelihood trace ----------
    hLL = figure(117); clf;
    plot(double(chain.ll(1:k)), 'k'); grid on;
    xlabel('Iteration'); ylabel('ll_{amp}+ll_{phase}');
    title(sprintf('Log-likelihood trace (k=%d)', k));
    set(gca,'FontSize',14);
    drawnow;

    saveas(hLL, sprintf('llcurrent_%d.jpg', k));
    % optional: savefig(hLL, sprintf('llcurrent_%d.fig', k), 'compact');
    close(hLL);

    % ---------- 2) parameter histograms (3x3) ----------
    hH = figure(114); clf;

    % (1) Temperature (C)
    subplot(3,3,1)
    histogram(mout(:,1)-273.15,'Normalization','Probability','BinWidth',10);
    xlabel('\circC'); xlim([700 1200]); title('Temperature');
    set(gca,'YTickLabel',[]);

    % (2) Gas flow rate
    subplot(3,3,2)
    histogram(mout(:,2),'Normalization','Probability','BinWidth',5);
    xlim([0 max(mout(:,2))+10]); xlabel('Kg/s'); title('Gas flow rate');
    set(gca,'YTickLabel',[]);

    % (3) Conduit radius
    subplot(3,3,3)
    histogram(mout(:,3),'Normalization','Probability','BinWidth',1);
    xlabel('m'); xlim([0 40]); title('Conduit Radius');
    set(gca,'YTickLabel',[]);

    % (4) Cap thickness
    subplot(3,3,4)
    histogram(mout(:,4),'Normalization','Probability','BinWidth',1);
    xlim([0 100]); xlabel('m'); title('Cap Thickness');
    set(gca,'YTickLabel',[]);

    % (5) Gas pocket thickness D = 10^(-m(6))
    subplot(3,3,5)
    D = 10.^(-mout(:,6));
    [~,edges] = histcounts(log10(D),50);
    histogram(D,10.^edges,'Normalization','Probability');
    xlabel('m'); xlim([1e-2 3]); title('Gas Pocket Thick');
    set(gca,'xscale','log'); xticks([1e-2 1e-1 1]); set(gca,'YTickLabel',[]);

    % (6) Permeability kappa = 10^(-m(5))
    subplot(3,3,6)
    kap = 10.^(-mout(:,5));
    [~,edges] = histcounts(log10(kap),50);
    histogram(kap,10.^edges,'Normalization','Probability');
    xlim([1e-12 5e-7]); xticks([1e-11 1e-9 1e-7]);
    xlabel('m^2'); title('Permeability');
    set(gca,'xscale','log'); set(gca,'YTickLabel',[]);

    % (7) Porosity phi(%) = 10^(m(7)) in your model convention
    subplot(3,3,7)
    por = 10.^(mout(:,7));
    [~,edges] = histcounts(log10(por),50);
    histogram(por,10.^edges,'Normalization','Probability');
    xlabel('%'); title('Porosity');
    set(gca,'xscale','log'); xlim([0.1 20]); xticks([0.2 2 20]);
    set(gca,'YTickLabel',[]);

    % (8) Density
    subplot(3,3,8)
    histogram(mout(:,8),'Normalization','Probability','BinWidth',3);
    xlabel('Kg/m^3'); title('Density');
    xlim([2800 3100]); set(gca,'YTickLabel',[]);

    % (9) Quality factor Qf
    subplot(3,3,9)
    histogram(mout(:,9),'Normalization','Probability','BinWidth',1);
    xlim([20 400]); title('Quality factor');
    set(gca,'YTickLabel',[]);

    drawnow;
    saveas(hH, sprintf('MCMC_source_%d.jpg', k));
    % optional: savefig(hH, sprintf('MCMC_source_%d.fig', k), 'compact');
    close(hH);

    % ---------- 3) spectral fit (MAP vs data) ----------
    % Use MAP.model which you already computed/stored in the loop.
    % Plot only inside your misfit band if desired; here I cap at 10 Hz like your old.
    hSpec = figure(118); clf;

    fmax = min(cfg.fPlotMax, 10);
    idxPlot = find(f <= fmax);

    plot(f(idxPlot), MAP.model.amp(idxPlot), 'r', ...
         f(idxPlot), LPspec.amp(idxPlot), 'b', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)','FontSize',18);
    ylabel('Amplitude','FontSize',18);
    title(sprintf('Ground displacement spectrum (MAP)  k=%d', k));
    xlim([0 fmax]);
    set(gca,'FontSize',18,'fontweight','bold');
    grid on;
    drawnow;

    saveas(hSpec, sprintf('fit_synthetic_%d.jpg', k));
    % optional: savefig(hSpec, sprintf('fit_synthetic_%d.fig', k), 'compact');
    close(hSpec);

    % ---------- 4) time-domain fit (MAP synthetic vs observed LP window) ----------
    % MAP.model.synthetic is the time series produced by forward model.
    % Apply the SAME filter as preprocessing for apples-to-apples.
    order=4; sampling_frequency=cfg.Fs; low_frequency=cfg.fLo; high_frequency=cfg.fHi;
    [b,a] = butter(order/2,[low_frequency high_frequency]/(sampling_frequency/2),'bandpass');

    syn = MAP.model.synthetic(:);
    syn_f = filtfilt(b,a,syn);

    obs = LP.x(:);
    obs_f = filtfilt(b,a,obs);

    dt = cfg.dt;
    t_syn = (0:numel(syn_f)-1)'*dt;
    t_obs = (0:numel(obs_f)-1)'*dt;

    % If lengths differ, align by truncation to common length
    Nmin = min(numel(syn_f), numel(obs_f));
    syn_f = syn_f(1:Nmin); obs_f = obs_f(1:Nmin);
    t_obs = t_obs(1:Nmin); t_syn = t_syn(1:Nmin);

    hTD = figure(122); clf;
    plot(t_syn, syn_f, 'r', t_obs, obs_f, 'b', 'LineWidth', 1.2);
    xlabel('time (s)','FontSize',18);
    ylabel('Amplitude (m)','FontSize',18);
    title(sprintf('Time-domain displacement (filtered)  k=%d', k));
    xlim([0 tau]);
    set(gca,'FontSize',18,'fontweight','bold');
    grid on;
    drawnow;

    saveas(hTD, sprintf('fit_timedomain_%d.jpg', k));
    % optional: savefig(hTD, sprintf('fit_timedomain_%d.fig', k), 'compact');
    close(hTD);
end