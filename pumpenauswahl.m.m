%% Pumpenauswahl: Versuchsstand Becken - Pumpe - Regelventil - Becken (Wasser)
% Platzhalterwerte in Abschnitt 2 und 3 durch eigene Daten ersetzen. Dann Run.
clear; clc; close all;

%% 1) Medium (Wasser, 20 °C)
rho  = 998;        % kg/m^3
nu   = 1.004e-6;   % m^2/s
pv   = 2339;       % Pa    Dampfdruck
patm = 101325;     % Pa
g    = 9.81;       % m/s^2

%% 2) Anlage (PLATZHALTER)
Q_soll    = 8;          % m^3/h  Soll-Volumenstrom
D         = 0.020;      % m      Rohr-Innendurchmesser
L         = 12;         % m      Gesamtrohrlaenge
eps_r     = 0.015e-3;   % m      Rauheit (Edelstahl)
zeta_E    = 6;          % -      Summe Einzelwiderstaende (Boegen, Filter, T)
zeta_V    = 4;          % -      Regelventil, aktuelle Stellung
zeta_Vmax = 15;         % -      max. sinnvolle Drosselung
H_geo     = 0;          % m      Hoehendifferenz (Kreislauf: 0)
h_s       = 0.3;        % m      Wasserstand ueber Pumpeneintritt (fuer NPSH)

%% 3) Pumpe: Datenblattpunkte (PLATZHALTER)
Qp    = [0    2     4     6     8     10    12];     % m^3/h
Hp    = [12.5 12.1  11.2  9.9   8.1   5.9   3.4];   % m
etap  = [0    0.22  0.38  0.47  0.50  0.45  0.35];  % -
NPSHp = [0.6  0.6   0.7   0.8   1.0   1.3   1.7];   % m

%% 4) Kennlinien auf Gitter
Q     = linspace(0.1, max(Qp), 400);                % m^3/h
HP    = interp1(Qp, Hp,    Q, 'pchip');
etaP  = interp1(Qp, etap,  Q, 'pchip');
NPSHr = interp1(Qp, NPSHp, Q, 'pchip');
HL    = Hloss(Q, D, L, eps_r, zeta_E, zeta_V, nu, g);
HA    = H_geo + HL;                                  % Anlagenkennlinie
NPSHa = (patm - pv)/(rho*g) + h_s - 0.5*HL;          % NPSH verfuegbar (halbe Verluste saugseitig angenommen)

%% 5) Betriebspunkt (Schnittpunkt HP = HA)
d = HP - HA;
i = find(d(1:end-1).*d(2:end) <= 0, 1, 'first');
if isempty(i)
    Q_op = NaN; H_op = NaN; eta_op = NaN; NPSHr_op = NaN; NPSHa_op = NaN;
else
    Q_op     = interp1(d(i:i+1), Q(i:i+1), 0);
    H_op     = interp1(Q, HP, Q_op);
    eta_op   = interp1(Q, etaP, Q_op);
    NPSHr_op = interp1(Q, NPSHr, Q_op);
    NPSHa_op = interp1(Q, NPSHa, Q_op);
end
P_el = rho*g*(Q_op/3600)*H_op/eta_op;                % W, elektrische Leistung

%% 6) Drosselstellung, die Q_soll einstellt
A        = pi*D^2/4;
vS       = (Q_soll/3600)/A;
HPs      = interp1(Qp, Hp, Q_soll, 'pchip');         % Pumpenhoehe bei Q_soll
Hf0      = Hloss(Q_soll, D, L, eps_r, zeta_E, 0, nu, g);
zeta_req = (HPs - H_geo - Hf0)/(vS^2/(2*g));
zV_plot  = max(zeta_req, 0);
HA_req   = H_geo + Hloss(Q, D, L, eps_r, zeta_E, zV_plot, nu, g);

%% 7) Urteil
ok1 = zeta_req >= 0;                 % Pumpe liefert genug Foerderhoehe bei Q_soll
ok2 = zeta_req <= zeta_Vmax;         % Pumpe nicht stark ueberdimensioniert
ok3 = NPSHa_op >= 1.3*NPSHr_op;      % NPSH-Reserve
ok4 = eta_op >= 0.7*max(etap);       % Wirkungsgrad im guten Bereich
if all([ok1 ok2 ok3 ok4])
    verdict = 'PUMPE PASST';
else
    verdict = 'PUMPE PASST NICHT';
end

fprintf('\n=== ERGEBNIS ===\n');
fprintf('Betriebspunkt (Ventil aktuell): Q = %.2f m^3/h, H = %.2f m\n', Q_op, H_op);
fprintf('Pumpe bei Q_soll: H = %.2f m\n', HPs);
fprintf('Noetige Drosselung: zeta_V = %.2f (max. %.1f)\n', zeta_req, zeta_Vmax);
fprintf('NPSH_A = %.2f m, NPSH_R = %.2f m, Reserve = %.2f (Ziel >= 1.3)\n', NPSHa_op, NPSHr_op, NPSHa_op/NPSHr_op);
fprintf('Wirkungsgrad = %.1f %%, P_el = %.0f W\n', 100*eta_op, P_el);
fprintf('Urteil: %s\n', verdict);

%% 8) Plot
figure('Color','w','Position',[100 100 1100 450]);

subplot(1,2,1); hold on; grid on; box on;
plot(Q, HP,    'b-',  'LineWidth', 2);
plot(Q, HA,    'r--', 'LineWidth', 1.5);
plot(Q, HA_req, 'm:', 'LineWidth', 1.8);
plot(Qp, Hp, 'bo', 'MarkerFaceColor', 'b');
plot(Q_op, H_op, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
plot([Q_soll Q_soll], [0 max(HP)], 'g-', 'LineWidth', 1.5);
xlabel('Q [m^3/h]'); ylabel('H [m]');
legend('Pumpe H_P(Q)', 'Anlage (Ventil aktuell)', 'Anlage (Drossel auf Q_{soll})', ...
       'Datenblatt', 'Betriebspunkt', 'Q_{soll}', 'Location', 'northeast');
title(sprintf('Urteil: %s', verdict));

subplot(1,2,2); hold on; grid on; box on;
yyaxis left
plot(Q, 100*etaP, 'g-', 'LineWidth', 1.8);
ylabel('\eta [%]');
yyaxis right
plot(Q, NPSHr, 'k-',  'LineWidth', 1.5);
plot(Q, NPSHa, 'k--', 'LineWidth', 1.5);
ylabel('NPSH [m]');
xlabel('Q [m^3/h]');
legend('\eta', 'NPSH_R', 'NPSH_A', 'Location', 'best');
title('Wirkungsgrad & NPSH');

%% Lokale Funktion: Verlusthoehe nach Darcy-Weisbach (Swamee-Jain, laminar <2300)
function H = Hloss(Qm3h, D, L, eps_r, zE, zV, nu, g)
    Qs  = Qm3h/3600;
    v   = Qs/(pi*D^2/4);
    Re  = max(v*D/nu, 1e-6);
    f   = 0.25./log10(eps_r/(3.7*D) + 5.74./Re.^0.9).^2;
    lam = Re < 2300;
    f(lam) = 64./Re(lam);
    H   = (f*L/D + zE + zV).*v.^2/(2*g);
end
