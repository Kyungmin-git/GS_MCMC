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