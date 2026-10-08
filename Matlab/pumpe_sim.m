% =========================================================================
% Pumpen-Simulation (MATLAB)
% Kennlinien, Arbeitspunkte, Wirkungsgrad, Anlauf (Euler), Temperaturregelung
% Aufruf: pumpe_sim  -> liest ../Pumpe_Daten.csv, schreibt Plots/ und Ergebnisse.csv
% Benoetigt: MATLAB R2020a+ (exportgraphics)
% =========================================================================
clear; close all; clc;

here = fileparts(mfilename('fullpath'));
plotDir = fullfile(here, 'Plots');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end

% ---------- Parameter aus CSV ----------
T = readtable(fullfile(here, '..', 'Pumpe_Daten.csv'), 'TextType', 'char');
W = str2double(string(T.Wert));           % Textwerte werden NaN und uebersprungen
P = struct();
for i = 1:height(T)
    if ~isnan(W(i))
        P.(matlab.lang.makeValidName(T.Parameter{i})) = W(i);
    end
end
A = pi * (P.D_innen / 2)^2;               % Rohrquerschnitt [m2]

% ---------- Modell (Funktionen als Handles) ----------
Hp    = @(Q) P.H0 * (1 - (Q ./ P.Q0).^2);                                  % Pumpenkennlinie [m]
eta   = @(Q) P.eta_max * exp(-((Q - P.Q_bep) ./ (0.5 * P.Q_bep)).^2);      % Wirkungsgrad [-]
Phyd  = @(Q, H) P.rho * P.g .* Q .* H;                                     % hydraulische Leistung [W]
Pel   = @(Q) Phyd(Q, Hp(Q)) ./ eta(Q);                                     % elektrische Leistung [W]
kges  = @(th) P.k_rohr + P.k_ventil0 .* (1 ./ th.^2 - 1);                  % Anlagenbeiwert [s2/m5]
Ha    = @(Q, th) P.H_geo + kges(th) .* Q.^2;                               % Anlagenkennlinie [m]
Qop   = @(th) sqrt((P.H0 - P.H_geo) ./ (P.H0 ./ P.Q0^2 + kges(th)));       % Arbeitspunkt Q [m3/s]

% ---------- Arbeitspunkt (Nennbetrieb, theta = 1) ----------
Q1  = Qop(1);
H1  = Hp(Q1);
Phy1 = Phyd(Q1, H1);
Pel1 = Pel(Q1);
eta1 = eta(Q1);

% ---------- Anlauf aus Ruhe ----------
dtA = 1e-4;
tA  = (0:round(0.5 / dtA))' * dtA;
QA  = zeros(size(tA));
kA  = P.g * A / P.L_rohr;
for i = 1:numel(tA) - 1
    QA(i+1) = QA(i) + dtA * kA * (Hp(QA(i)) - Ha(QA(i), 1));
end
idx95 = find(QA >= 0.95 * Q1, 1, 'first');
t95   = tA(idx95);

% ---------- Beckentemperatur (Hysterese, 1 h, dt = 1 s) ----------
Qverl = Pel(Q1) - Phy1;                  % Pumpenverlust ins Wasser [W]
C     = P.rho * P.V_tank * P.c_w;        % Waermekapazitaet [J/K]
ts    = P.T_soll;
hy    = P.Hysterese;
dtT   = 1;
tT    = (0:round(3600 / dtT))' * dtT;
nT    = numel(tT);
TT    = zeros(nT, 1);
TT(1) = P.T_start;
heiz  = false(nT, 1);
kuehl = false(nT, 1);
han = false; kan = false;
for i = 1:nT - 1
    Ti = TT(i);
    if Ti < ts - hy
        han = true;
    elseif Ti >= ts
        han = false;
    end
    if Ti > ts + hy
        kan = true;
    elseif Ti <= ts
        kan = false;
    end
    q = han * P.P_heiz - kan * P.P_kuehl + Qverl - P.UA * (Ti - P.T_umg);
    TT(i+1) = Ti + dtT * q / C;
    heiz(i) = han;
    kuehl(i) = kan;
end
heiz(nT) = han;
kuehl(nT) = kan;

% ---------- Ergebnisse ----------
names = ["Foerderstrom_Q_op"; "Foerderhoehe_H_op"; "Hydraulische_Leistung"; ...
          "Wirkungsgrad_op"; "Elektrische_Leistung_op"; "Verlustleistung_Pumpe_ins_Wasser"; ...
          "Anlaufzeit_t95"; "Heizelement_Anteil"; "Kuehlung_Anteil"; ...
          "Temperatur_Ende"; "Temperatur_Max"];
werte = [Q1*6e4; H1; Phy1; eta1; Pel1; Qverl; t95; mean(heiz)*100; mean(kuehl)*100; ...
         TT(end); max(TT)];
einh  = ["L/min"; "m"; "W"; "-"; "W"; "W"; "s"; "%"; "%"; "C"; "C"];
ERG   = table(names, werte, einh, 'VariableNames', {'Groesse', 'Wert', 'Einheit'});
writetable(ERG, fullfile(here, 'Ergebnisse.csv'));
disp(ERG);

% ---------- Plots ----------
Qs = linspace(0, P.Q0, 400);
thetas = [1.0 0.75 0.5 0.3];

% 1 Kennlinien + Arbeitspunkte
f = figure('Color', 'w', 'Position', [100 100 800 500]);
hold on; grid on; box on;
plot(Qs * 6e4, Hp(Qs), 'k-', 'LineWidth', 2, 'DisplayName', 'Pumpenkennlinie');
for th = thetas
    plot(Qs * 6e4, Ha(Qs, th), '--', 'LineWidth', 1.2, 'DisplayName', sprintf('Anlage theta=%.2f', th));
    Qo = Qop(th);
    plot(Qo * 6e4, Hp(Qo), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'auto', ...
        'DisplayName', sprintf('AP theta=%.2f: %.1f L/min, %.1f m', th, Qo*6e4, Hp(Qo)));
end
xlabel('Foerderstrom Q [L/min]'); ylabel('Foerderhoehe H [m]');
title('Pumpen- und Anlagenkennlinie, Arbeitspunkte');
legend('Location', 'northeast');
exportgraphics(f, fullfile(plotDir, '1_Kennlinien.png'), 'Resolution', 150);
close(f);

% 2 Wirkungsgrad + elektrische Leistung
f = figure('Color', 'w', 'Position', [100 100 800 700]);
subplot(2, 1, 1); hold on; grid on; box on;
plot(Qs * 6e4, eta(Qs), 'g-', 'LineWidth', 2);
xline(P.Q_bep * 6e4, ':', 'Bestpunkt');
ylabel('Wirkungsgrad eta [-]');
subplot(2, 1, 2); hold on; grid on; box on;
plot(Qs * 6e4, Pel(Qs), 'r-', 'LineWidth', 2, 'DisplayName', 'P_{el}');
yline(P.P_nenn, 'k--', 'P_{nenn}', 'DisplayName', 'P_nenn');
xlabel('Foerderstrom Q [L/min]'); ylabel('Elektrische Leistung [W]');
legend('Location', 'northwest');
sgtitle('Wirkungsgrad und Leistungsaufnahme');
exportgraphics(f, fullfile(plotDir, '2_Wirkungsgrad_Leistung.png'), 'Resolution', 150);
close(f);

% 3 Ventilstellung -> Arbeitspunkt
th = linspace(0.2, 1.0, 60);
Qv = Qop(th);
Hv = Hp(Qv);
f = figure('Color', 'w', 'Position', [100 100 800 900]);
subplot(3, 1, 1); plot(th * 100, Qv * 6e4, 'b-', 'LineWidth', 2); grid on; box on; ylabel('Q [L/min]');
subplot(3, 1, 2); plot(th * 100, Hv, 'm-', 'LineWidth', 2); grid on; box on; ylabel('H [m]');
subplot(3, 1, 3); plot(th * 100, Pel(Qv), 'r-', 'LineWidth', 2); grid on; box on;
ylabel('P_{el} [W]'); xlabel('Ventilstellung theta [%]');
sgtitle('Drosselventil: Arbeitspunkt in Abhaengigkeit der Stellung');
exportgraphics(f, fullfile(plotDir, '3_Ventilstellung.png'), 'Resolution', 150);
close(f);

% 4 Anlauf
f = figure('Color', 'w', 'Position', [100 100 800 500]);
hold on; grid on; box on;
plot(tA * 1000, QA * 6e4, 'b-', 'LineWidth', 2, 'DisplayName', 'Q(t)');
yline(Q1 * 6e4, 'k--', 'DisplayName', 'Q_{op}');
xline(t95 * 1000, ':', sprintf('t95 = %.0f ms', t95*1000), 'DisplayName', 't95');
xlabel('Zeit t [ms]'); ylabel('Foerderstrom Q [L/min]');
title('Anlauf aus Ruhe (Euler, dt = 0.1 ms)');
legend('Location', 'southeast');
exportgraphics(f, fullfile(plotDir, '4_Anlauf.png'), 'Resolution', 150);
close(f);

% 5 Temperatur
f = figure('Color', 'w', 'Position', [100 100 800 700]);
subplot(2, 1, 1); hold on; grid on; box on;
patch([0 60 60 0], [ts-hy ts-hy ts+hy ts+hy], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(tT / 60, TT, 'r-', 'LineWidth', 2, 'DisplayName', 'T Becken');
yline(ts, 'k--', 'DisplayName', 'T_{soll}');
ylabel('Temperatur [C]');
legend('Location', 'southeast');
subplot(2, 1, 2); hold on; grid on; box on;
stairs(tT / 60, double(heiz), 'Color', [1 0.6 0], 'LineWidth', 2, 'DisplayName', 'Heizen');
stairs(tT / 60, -double(kuehl), 'b-', 'LineWidth', 2, 'DisplayName', 'Kuehlen');
xlabel('Zeit t [min]'); ylabel('Zustand'); yticks([-1 0 1]);
legend('Location', 'southwest');
sgtitle('Beckentemperatur mit Hysterese-Regler (1 h)');
exportgraphics(f, fullfile(plotDir, '5_Temperatur.png'), 'Resolution', 150);
close(f);

fprintf('Plots gespeichert in: %s\n', plotDir);
