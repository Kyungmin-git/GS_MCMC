
function S = compute_spectrum(x, dt, Fs)
    % FFT normalization: keep consistent with your forward model convention.
    % Here we keep a simple, documented convention:
    %   X = fft(x)/sqrt(N) * sqrt(dt/df)
    N = numel(x);
    df = (1/dt)/N;

    X = fft(x(:), N) / sqrt(N) * sqrt(dt/df);
    f = (0:(N/2))' * (Fs/N);

    Xpos = X(1:N/2+1);

    S.f = f;
    S.df = df;
    S.complex = Xpos;
    S.amp = abs(Xpos);
    S.phase = angle(Xpos);
end