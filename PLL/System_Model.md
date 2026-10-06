%% GNSS PLL PARAMETER SWEEP  (GPS L1 / L2 / L5)
% Sweeps a charge-pump PLL synthesizer over loop bandwidth, phase margin,
% charge-pump current, reference divider (step frequency) and integer- vs
% fractional-N, then ranks the designs for each GNSS band on:
%   RMS phase jitter, reference spur, overshoot, lock time, cap area.
%
% Based on the phase-domain method of the MathWorks example
% "Model PLLs in the Phase Domain", but self-contained:
%   - needs only the Control System Toolbox
%   - no getNoiseTransferImpedance / getFlickerFilter helpers needed
%   - the loop filter (3rd-order passive, same topology as the example)
%     is DESIGNED from the target bandwidth + phase margin for every point,
%     instead of using fixed R/C values.
%
% Topology (charge pump injects at node 1):
%   node1 --C1-- gnd
%   node1 --R2--C2-- gnd           (stabilizing zero)
%   node1 --R3-- node3 --C3-- gnd  (spur-suppression pole), Vtune = V(node3)
%
% Outputs: results table (workspace + CSV), best design per band,
%          plots of phase noise contributions, step response, trade-offs.

clear; clc; close all;

%% ===================== 1. GNSS BANDS ====================================
% sigBW   : two-sided signal bandwidth the IF must accommodate
% IFtarget: preferred IF (script picks the nearest feasible one)
bands = struct( ...
  'name',    {'L1',        'L2',       'L5'}, ...
  'fRF',     {1575.42e6,   1227.60e6,  1176.45e6}, ...
  'sigBW',   {2.046e6,     2.046e6,    20.46e6}, ...  % C/A, L2C, L5 main lobe
  'IFtarget',{4.092e6,     4.092e6,    14.322e6});
% Optional extra bands (uncomment to add):
% bands(end+1) = struct('name','L3','fRF',1381.05e6, 'sigBW',2.046e6,'IFtarget',4.092e6);
% bands(end+1) = struct('name','L4','fRF',1379.913e6,'sigBW',2.046e6,'IFtarget',4.092e6);

IFmax      = 20e6;   % highest IF your ADC/IF chain supports (Hz)
allowHighSideLO = false;   % true -> also consider LO above RF (negative IF)

%% ===================== 2. FIXED HARDWARE ASSUMPTIONS ====================
P.fRef    = 16.368e6;  % TCXO reference (Hz)
P.Kvco    = 40e6;      % VCO gain (Hz/V) - replace with your VCO's value
P.k3      = 6;         % 3rd pole placed at k3 * loop bandwidth
P.kC3     = 10;        % C3 = C1 / kC3
P.modBits = 24;        % fractional-N modulus = 2^modBits
P.mashOrd = 3;         % MASH sigma-delta order (fractional-N only)
P.Temp    = 298;       % K

% Noise models -----------------------------------------------------------
% TCXO phase noise (dBc/Hz vs offset Hz) - typical 16 MHz TCXO
P.refPN.f = [10   100  1e3  1e4  1e5  1e6];
P.refPN.L = [-90 -120 -140 -150 -155 -155];
% Free-running VCO phase noise, specified at P.vcoPNfreq; scaled
% 20*log10(fLO/P.vcoPNfreq) for other bands
P.vcoPN.f = [1e4  1e5  1e6  1e7];
P.vcoPN.L = [-82 -107 -128 -148];
P.vcoPNfreq = 1.5713e9;
% PFD/CP/divider in-band noise via figure of merit
P.FOM      = -222;   % normalized flat FOM (dBc/Hz), output = FOM+10log(fPFD)+20log(N)
P.FOMflick = -120;   % normalized 1/f noise at 10 kHz offset, 1 GHz carrier (dBc/Hz)

% Reference-spur model ---------------------------------------------------
P.Ileak    = 1e-9;   % leakage at Vtune node (A) - raise this for FeCAP leakage
P.cpMismatch = 0.02; % UP/DN current mismatch (fraction of Icp)
P.tReset   = 200e-12;% PFD reset (anti-backlash) pulse width (s)

% Lock-time test ----------------------------------------------------------
P.fStep  = 1e6;      % output frequency step for lock-time test (Hz)
P.fTol   = 100;      % "locked" when within this error of final (Hz)

% Jitter integration -----------------------------------------------------
P.fIntLow = 10;      % ~carrier tracking loop bandwidth (Hz)
% upper limit = sigBW/2 of each band (capped at fPFD/2)

%% ===================== 3. SWEEP RANGES ==================================
sweep.fc    = [30 50 100 150 200 300 500 800]*1e3;  % loop bandwidth (Hz)
sweep.PM    = [45 55 65];                           % phase margin (deg)
sweep.Icp   = [0.1 0.25 0.5 1 2]*1e-3;              % charge pump (A)
sweep.Rdiv  = [1 2 4];                              % reference divider
sweep.mode  = {'integer','fractional'};

%% ===================== 4. PASS/FAIL LIMITS AND RANKING ==================
lim.PMmin      = 40;     % deg (measured; lands ~1 deg below the target)
lim.jitterMax  = 3;      % deg RMS
lim.spurMax    = -70;    % dBc
lim.OSmax      = 30;     % % overshoot
lim.lockMax    = 100e-6; % s
lim.CtotMax    = 500e-12;% F total loop-filter capacitance (chip area!)
% Cost weights among passing designs (each metric normalized to its limit)
wgt.jitter = 3; wgt.spur = 1; wgt.OS = 1; wgt.lock = 1; wgt.Ctot = 2;

%% ===================== 5. RUN THE SWEEP =================================
[gB,gM,gR,gI,gF,gP] = ndgrid(1:numel(bands), 1:numel(sweep.mode), ...
    1:numel(sweep.Rdiv), 1:numel(sweep.Icp), 1:numel(sweep.fc), 1:numel(sweep.PM));
nRun = numel(gB);
rows = cell(nRun,1);
fprintf('Running %d designs...\n', nRun);
tic;
for k = 1:nRun
    b    = bands(gB(k));
    mode = sweep.mode{gM(k)};
    rows{k} = evaluatePLL(b, mode, sweep.Rdiv(gR(k)), sweep.Icp(gI(k)), ...
                          sweep.fc(gF(k)), sweep.PM(gP(k)), P, IFmax, allowHighSideLO, false);
    if mod(k, ceil(nRun/10)) == 0
        fprintf('  %3.0f%%  (%.1f s)\n', 100*k/nRun, toc);
    end
end
rows = rows(~cellfun(@isempty, rows));
res  = struct2table([rows{:}]');

%% ===================== 6. PASS/FAIL AND COST ============================
res.Pass = res.PM_deg >= lim.PMmin & res.Jitter_deg <= lim.jitterMax & ...
           res.RefSpur_dBc <= lim.spurMax & res.Overshoot_pct <= lim.OSmax & ...
           res.LockTime_us <= lim.lockMax*1e6 & res.Ctot_pF <= lim.CtotMax*1e12 & ...
           res.Linear & res.BWvalid;
res.Cost = wgt.jitter*res.Jitter_deg/lim.jitterMax + ...
           wgt.spur  *10.^((res.RefSpur_dBc - lim.spurMax)/20) + ...
           wgt.OS    *res.Overshoot_pct/lim.OSmax + ...
           wgt.lock  *res.LockTime_us/(lim.lockMax*1e6) + ...
           wgt.Ctot  *res.Ctot_pF/(lim.CtotMax*1e12);
res.Cost(~res.Pass) = Inf;
res = sortrows(res, {'Band','Cost'});

writetable(res, 'gnss_pll_sweep_results.csv');
save('gnss_pll_sweep_results.mat', 'res', 'P', 'sweep', 'lim', 'wgt', 'bands');
fprintf('\n%d designs evaluated, %d pass. Saved gnss_pll_sweep_results.csv/.mat\n', ...
        height(res), nnz(res.Pass));

%% ===================== 7. REPORT BEST DESIGN PER BAND ===================
showCols = {'Band','Mode','Rdiv','Step_Hz','N','IF_MHz','Icp_mA','fc_kHz','PM_deg', ...
            'Jitter_deg','RefSpur_dBc','Overshoot_pct','LockTime_us','Ctot_pF'};
for ib = 1:numel(bands)
    sub = res(strcmp(res.Band, bands(ib).name) & res.Pass, :);
    fprintf('\n===== %s (%.3f MHz): %d passing designs =====\n', ...
            bands(ib).name, bands(ib).fRF/1e6, height(sub));
    if isempty(sub)
        fprintf('No design meets all limits - relax lim.* or widen the sweep.\n');
        continue
    end
    disp(sub(1:min(5,height(sub)), showCols));
    best = sub(1,:);
    fprintf('Best: R2=%.3g ohm, R3=%.3g ohm, C1=%.3g pF, C2=%.3g pF, C3=%.3g pF\n', ...
        best.R2, best.R3, best.C1_pF, best.C2_pF, best.C3_pF);

    % Re-run the winner with full detail for plotting
    d = evaluatePLL(bands(ib), best.Mode{1}, best.Rdiv, best.Icp_mA*1e-3, ...
                    best.fcTarget_kHz*1e3, best.PMtarget, P, IFmax, allowHighSideLO, true);
    plotBest(d, bands(ib).name);
end

%% ===================== 8. TRADE-OFF PLOTS ==============================
for ib = 1:numel(bands)
    sub = res(strcmp(res.Band, bands(ib).name) & strcmp(res.Mode,'integer') & ...
              res.PMtarget == 55 & res.Rdiv == min(sweep.Rdiv), :);
    if isempty(sub), continue; end
    figure('Name', [bands(ib).name ' trade-offs']);
    tiledlayout(1,3);
    for m = 1:3
        nexttile; hold on; grid on;
        for ii = 1:numel(sweep.Icp)
            s2 = sub(abs(sub.Icp_mA - sweep.Icp(ii)*1e3) < 1e-9, :);
            switch m
                case 1, y = s2.Jitter_deg;  yl = 'RMS phase jitter (deg)';
                case 2, y = s2.RefSpur_dBc; yl = 'Reference spur (dBc)';
                case 3, y = s2.Ctot_pF;     yl = 'Total loop-filter C (pF)';
            end
            semilogx(s2.fc_kHz, y, '-o', 'DisplayName', sprintf('I_{CP}=%g mA', sweep.Icp(ii)*1e3));
        end
        set(gca,'XScale','log'); xlabel('Measured loop bandwidth (kHz)'); ylabel(yl);
        if m == 3, set(gca,'YScale','log'); legend('Location','best'); end
    end
    sgtitle(sprintf('%s, integer-N, PM target 55%c, R=%d', bands(ib).name, char(176), min(sweep.Rdiv)));
end

%% =======================================================================
%%                           LOCAL FUNCTIONS
%% =======================================================================

function out = evaluatePLL(b, mode, Rdiv, Icp, fcT, PMt, P, IFmax, allowHigh, detail)
% Designs and evaluates one PLL. Returns [] if no valid frequency plan.
out = [];
kB   = 1.380649e-23;
fPFD = P.fRef / Rdiv;

% ---- Frequency plan / step frequency ----
switch mode
    case 'integer'
        [N, IF] = pickIntegerN(b, fPFD, IFmax, allowHigh);
        if isempty(N), return; end
        stepHz = fPFD;
        bndOff = NaN;
    case 'fractional'
        IF     = b.IFtarget;
        N      = (b.fRF - IF) / fPFD;
        stepHz = fPFD / 2^P.modBits;
        fr     = N - floor(N);
        bndOff = min(fr, 1-fr) * fPFD;  % integer-boundary spur offset (Hz)
end
fLO = b.fRF - IF;

% ---- Loop filter design (2nd order + extra pole) ----
fp3 = P.k3 * fcT;
wc  = 2*pi*fcT;
phi = deg2rad(PMt + atand(fcT/fp3));   % pre-compensate the phase the 3rd pole eats
T1  = (sec(phi) - tan(phi)) / wc;
T2  = 1 / (wc^2 * T1);
A0  = Icp*P.Kvco/(N*wc^2) * sqrt((1+(wc*T2)^2)/(1+(wc*T1)^2));
C1  = A0*T1/T2;  C2 = A0 - C1;  R2 = T2/C2;
C3  = C1 / P.kC3; R3 = 1/(2*pi*fp3*C3);

% ---- Models ----
Zlf = loopFilterSS(R2, R3, C1, C2, C3);            % inputs [Icp, iR2, iR3] -> Vtune
vco = ss(0, 2*pi*P.Kvco, 1, 0);
G   = (Icp/(2*pi)) * Zlf(1,1) * vco;               % forward gain
Lg  = G / N;                                       % loop gain
[~, PMm, ~, wcp] = margin(Lg);
fcM = wcp/(2*pi);

% ---- Phase noise over offset frequency ----
fHi = min(b.sigBW/2, fPFD/2);
f   = logspace(0, log10(fPFD/2), 1500);  w = 2*pi*f;
Zf  = squeeze(freqresp(Zlf, w));                   % 3 x nf
Kv  = 2*pi*P.Kvco ./ (1j*w);
Lf  = (Icp/(2*pi)) * Zf(1,:) .* Kv / N;
Hf  = N * Lf ./ (1 + Lf);                          % ref phase -> out phase
Ef  = 1 ./ (1 + Lf);                               % VCO phase -> out phase

S.ref = pn2S(P.refPN.f, P.refPN.L, f) .* abs(Hf).^2;
S.vco = pn2S(P.vcoPN.f, P.vcoPN.L + 20*log10(fLO/P.vcoPNfreq), f) .* abs(Ef).^2;
Sin   = 2*10^((P.FOM + 10*log10(fPFD))/10) + ...
        2*10.^((P.FOMflick + 20*log10(fLO/1e9) - 10*log10(f/1e4))/10) / N^2;
if strcmp(mode,'fractional')                       % MASH quantization noise
    Sin = Sin + (2*pi)^2/(12*fPFD) * (2*sin(pi*f/fPFD)).^(2*(P.mashOrd-1));
end
S.pll = Sin .* abs(Hf).^2;
S.R2  = 4*kB*P.Temp/R2 * abs(Zf(2,:)).^2 .* abs(Kv.*Ef).^2;
S.R3  = 4*kB*P.Temp/R3 * abs(Zf(3,:)).^2 .* abs(Kv.*Ef).^2;
S.tot = S.ref + S.vco + S.pll + S.R2 + S.R3;       % one-sided rad^2/Hz

m = f >= P.fIntLow & f <= fHi;
sigPhi = sqrt(trapz(f(m), S.tot(m)));              % rad RMS

% ---- Reference spur (narrowband FM from CP ripple at fPFD) ----
Zp   = freqresp(Zlf(1,1), 2*pi*fPFD);
Ih   = 2*(P.Ileak + P.cpMismatch*Icp*P.tReset*fPFD);
beta = P.Kvco * Ih * abs(Zp) / fPFD;
spur = 20*log10(beta/2);

% ---- Frequency-step response: overshoot and lock time ----
Tcl  = feedback(G, 1/N) / N;                       % DC gain 1
tEnd = 300/(2*pi*fcT);
t    = linspace(0, tEnd, 20000)';
y    = step(Tcl, t);
OS   = max(0, (max(y) - 1)*100);
idx  = find(abs(y - 1) > P.fTol/P.fStep, 1, 'last');
if isempty(idx), tLock = 0; elseif idx == numel(t), tLock = Inf; else, tLock = t(idx+1); end

% Peak PFD phase error during the step: linear model valid only if < 2*pi
Eint = ss(0,1,1,0) * feedback(ss(1), Lg);
pe   = step(Eint, t) * 2*pi*P.fStep/N;
peak = max(abs(pe));

% ---- Pack results ----
out = struct('Band', b.name, 'Mode', mode, 'Rdiv', Rdiv, 'fPFD_MHz', fPFD/1e6, ...
    'Step_Hz', stepHz, 'N', N, 'IF_MHz', IF/1e6, 'fLO_MHz', fLO/1e6, ...
    'BoundarySpurOffset_kHz', bndOff/1e3, 'Icp_mA', Icp*1e3, ...
    'fcTarget_kHz', fcT/1e3, 'PMtarget', PMt, 'fc_kHz', fcM/1e3, 'PM_deg', PMm, ...
    'Jitter_deg', rad2deg(sigPhi), 'Jitter_fs', sigPhi/(2*pi*fLO)*1e15, ...
    'RefSpur_dBc', spur, 'Overshoot_pct', OS, 'LockTime_us', tLock*1e6, ...
    'PeakPhaseErr_rad', peak, 'Linear', peak < 2*pi, 'BWvalid', fcM <= fPFD/10, ...
    'R2', R2, 'R3', R3, 'C1_pF', C1*1e12, 'C2_pF', C2*1e12, 'C3_pF', C3*1e12, ...
    'Ctot_pF', (C1+C2+C3)*1e12);
if detail
    out.f = f; out.S = S; out.t = t; out.y = y; out.fHi = fHi; out.fIntLow = P.fIntLow;
    out.fStep = P.fStep;
end
end

function Z = loopFilterSS(R2, R3, C1, C2, C3)
% State space of the 3rd-order passive filter.
% States [v1; vC2; v3], inputs [i_CP; i_noise_R2; i_noise_R3], output v3.
% Noise currents are Norton sources in parallel with each resistor.
% (To model FeCAP leakage later, add -v/Rleak terms to the diagonal of A.)
A = [-(1/R2 + 1/R3)/C1,  1/(R2*C1),  1/(R3*C1);
       1/(R2*C2),       -1/(R2*C2),  0;
       1/(R3*C3),        0,         -1/(R3*C3)];
B = [1/C1, -1/C1, -1/C1;
     0,     1/C2,  0;
     0,     0,     1/C3];
Z = ss(A, B, [0 0 1], zeros(1,3));
end

function [N, IF] = pickIntegerN(b, fPFD, IFmax, allowHigh)
% Integer N whose IF is wide enough for the signal and closest to IFtarget.
IFmin = b.sigBW/2;
Ncand = floor((b.fRF - IFmax)/fPFD) : ceil((b.fRF + IFmax)/fPFD);
IFc   = b.fRF - Ncand*fPFD;                 % >0 low-side LO, <0 high-side LO
ok    = abs(IFc) >= IFmin & abs(IFc) <= IFmax;
if ~allowHigh, ok = ok & IFc > 0; end
if ~any(ok), N = []; IF = []; return; end
Ncand = Ncand(ok); IFc = IFc(ok);
[~, i] = min(abs(abs(IFc) - b.IFtarget));
N = Ncand(i); IF = IFc(i);
end

function S = pn2S(fo, Ldb, f)
% SSB phase noise table (dBc/Hz) -> one-sided phase PSD (rad^2/Hz).
% Log-linear interpolation, slope-extrapolated below the table, flat above.
x  = interp1(log10(fo), Ldb, log10(f), 'linear');
lo = f < fo(1);
sl = (Ldb(2)-Ldb(1)) / (log10(fo(2))-log10(fo(1)));
x(lo) = Ldb(1) + sl*(log10(f(lo)) - log10(fo(1)));
x(f > fo(end)) = Ldb(end);
S = 2 * 10.^(x/10);
end

function plotBest(d, bandName)
figure('Name', [bandName ' best design']);
tiledlayout(1,2);
nexttile;
toL = @(S) 10*log10(S/2);
semilogx(d.f, toL(d.S.tot), 'k', 'LineWidth', 2); hold on;
semilogx(d.f, toL(d.S.ref)); semilogx(d.f, toL(d.S.vco)); semilogx(d.f, toL(d.S.pll));
semilogx(d.f, toL(d.S.R2));  semilogx(d.f, toL(d.S.R3));
xline(d.fIntLow, ':'); xline(d.fHi, ':'); xline(d.fc_kHz*1e3, '--', 'f_c');
grid on; ylim([-180 -60]);
xlabel('Offset (Hz)'); ylabel('L(f) (dBc/Hz)');
legend('Total','TCXO','VCO','PFD/CP/div','R2','R3','Location','southwest');
title(sprintf('%s %s-N: %.2f%c RMS, spur %.0f dBc', bandName, d.Mode, ...
      d.Jitter_deg, char(176), d.RefSpur_dBc));
nexttile;
plot(d.t*1e6, d.y*d.fStep/1e6); grid on;
xlabel('Time (\mus)'); ylabel('Output frequency change (MHz)');
title(sprintf('f_c=%.0f kHz, PM=%.1f%c, OS=%.1f%%, lock=%.1f \\mus', ...
      d.fc_kHz, d.PM_deg, char(176), d.Overshoot_pct, d.LockTime_us));
end
