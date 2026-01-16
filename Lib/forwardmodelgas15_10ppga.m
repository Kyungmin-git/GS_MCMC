function [Posc,synthetic,synangle,uspect,A_res,A_exc,A_path,uspectf,anglef] = forwardmodelgas15_10ppga(m,t0,qn,sigma_k_log)

%% ---------------- PHYSICAL PARAMETERS -----------------------------------
mu_g  = 1e-5;        % gas viscosity (Pa s)
M     = 0.018;       % molecular weight (kg/mol), water vapor
Rg    = 8.3145;      % J/(mol K)
Pex   = 101325;      % Pa (atmospheric)

T     = m(1);
Q     = m(2);
R     = m(3);
L     = m(4);
kappa = 10^(-m(5));
D     = 10^(-m(6));
phi   = 0.01*10^(m(7));   % convert from % to fraction
rho_s = m(8);
Qf    = m(9);

S = pi*R^2;

%% ---------------- AUXILIARY PARAMETERS ----------------------------------
den = (Pex - Rg*T*Q^2/(S^2*phi^2*M*Pex));  % common denominator (document in paper)
beta_a = S*phi*M/(Rg*T);
beta_b = mu_g*phi/(kappa*den);
beta_c = Pex*M/(Rg*T*den);
beta_d = 2*Q/(S*phi*den);
beta_e = S*M*D/(Rg*T);
P0     = Pex + mu_g*Rg*T*Q*L/(S*kappa*M*den);

%% ---------------- HARMONIC OSCILLATOR COEFFICIENTS ----------------------
GAMMA0 = 1;
GAMMA1 = (2*(beta_a*beta_d + beta_b*beta_e)*L + beta_a*beta_b*L^2)/(2*beta_a);
GAMMA2 = (2*beta_c*beta_e*L + beta_a*beta_c*L^2)/(2*beta_a);

gamma0 = beta_b*L/beta_a;
gamma1 = beta_c*L/beta_a;

%% ---------------- SIMULATION TIME GRID ----------------------------------
global LPstartt LPendt distance

tau = LPendt - LPstartt;

Fs_model = 50;             % Hz (matches your data)
dt = 1/Fs_model;           % s
time1 = 0:dt:tau;
Lw = numel(time1);
df = (1/dt)/Lw;

%% ---------------- SOURCE/PATH TRANSFER FUNCTIONS ------------------------
i = 1i;

u_z   = zeros(1,Lw);
A_p   = zeros(1,Lw);
A_res = zeros(1,Lw);
A_exc = zeros(1,Lw);
A_path= zeros(1,Lw);

N_exc = numel(qn);

sigma_k = 10^(sigma_k_log);
sigma_row = sigma_k * (qn(:)'.^(1/3));     % 1xN_exc



mm = 1;

for j = 0:(1/tau):(1/(2*dt))   % positive freqs (including 0)
    ome = 2*pi*j;

    % reservoir/source response
    A_res(mm) = (gamma0*(i*ome)^0 + gamma1*(i*ome)^1) / ...
                (GAMMA0*(i*ome)^0 + GAMMA1*(i*ome)^1 + GAMMA2*(i*ome)^2);

    % excitation: sum of Gaussian impulses
    A_exc(mm) = sum( qn .* exp(-i*ome*t0) .* exp(-0.5*(sigma_row.^2)*(ome^2)) );

    A_p(mm) = A_res(mm)*A_exc(mm);

    % path term
    if ome == 0
        A_path(mm) = 0;
    else
        f = ome/(2*pi);
        vc = 1295 * f^(-0.374);    % phase velocity
        vu = 0.73 * vc;            % group velocity
        A_path(mm) = exp(i*((ome/vc)*distance + pi/4)) ...
            * (1/(8*rho_s*vc^2*vu)) ...
            * sqrt(2*vc*ome/(pi*distance)) ...
            * exp(-ome*distance/(2*vu*Qf));
    end

    % vertical ground displacement spectrum
    u_z(mm) = S * A_res(mm) * A_exc(mm) * A_path(mm);

    mm = mm + 1;
end

%% ---------------- BACK TO TIME DOMAIN -----------------------------------
% Pressure evolution (steady-state frequency approach)
P = ifft(A_p,'symmetric') * sqrt(Lw)/sqrt(dt/df);
Posc = P + Pex - P0;

% Synthetic displacement
U_z = ifft(u_z,'symmetric') * sqrt(Lw)/sqrt(dt/df);
synthetic = real(U_z);

% Unfiltered spectrum outputs
uspect   = abs(u_z);
synangle = angle(u_z);

%% ---------------- FILTERED SPECTRUM (1–20 Hz) ---------------------------
order = 4;
low_frequency  = 1;
high_frequency = 20;

[b,a] = butter(order/2, [low_frequency high_frequency]/(Fs_model/2), 'bandpass');
syntheticf = filtfilt(b,a,synthetic);

% Recompute FFT scaling CONSISTENTLY for syntheticf length
Lwf = numel(syntheticf);
df_f = (1/dt)/Lwf;

U = fft(syntheticf(:), Lwf) / sqrt(Lwf) * sqrt(dt/df_f);

% positive frequencies only
Upos = U(1:(Lwf/2+1));

uspectf = abs(Upos);
anglef  = angle(Upos);

end
