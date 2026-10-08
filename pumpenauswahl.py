# -*- coding: utf-8 -*-
"""Pumpenauswahl: Versuchsstand Becken - Pumpe - Regelventil - Becken (Wasser)

Interaktiv im Terminal: Werte eingeben, Enter = Standardwert.
Listen (Pumpendaten) mit ';' trennen. Dezimal: '.' oder ','.
Benoetigt: numpy, matplotlib
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


# ---------- Eingabe ----------
def ask(text, default):
    """Zahl abfragen; Enter = Standardwert."""
    while True:
        s = input(f"{text} [{default}]: ").strip().replace(",", ".")
        if s == "":
            return float(default)
        try:
            return float(s)
        except ValueError:
            print("  Bitte eine Zahl eingeben.")


def ask_list(text, default):
    """Werteliste abfragen (Trennzeichen ';'); Enter = Standardwert."""
    default_str = "; ".join(str(x) for x in default)
    while True:
        s = input(f"{text} [{default_str}]: ").strip()
        if s == "":
            return np.array(default, dtype=float)
        try:
            arr = np.array([float(x.strip().replace(",", ".")) for x in s.split(";")])
        except ValueError:
            print("  Nur Zahlen, getrennt durch ';'.")
            continue
        if arr.size < 3:
            print("  Mindestens 3 Werte.")
            continue
        return arr


print("=== EINGABE (Enter = Standardwert) ===")
print("-- Medium: Wasser 20 °C (fest)")
print("-- Anlage")
Q_soll = ask("Soll-Volumenstrom Q_soll [m^3/h]", 8.0)
D = ask("Rohr-Innendurchmesser D [mm]", 20.0) / 1000.0
L = ask("Gesamtrohrlaenge L [m]", 12.0)
eps_r = ask("Rauheit eps [mm] (Edelstahl ca. 0.015)", 0.015) / 1000.0
zeta_E = ask("Summe Einzelwiderstaende zeta_E [-]", 6.0)
zeta_V = ask("Regelventil zeta_V, aktuelle Stellung [-]", 4.0)
zeta_Vmax = ask("max. sinnvolle Drosselung zeta_Vmax [-]", 15.0)
H_geo = ask("Hoehendifferenz Becken->Pumpe H_geo [m] (Kreislauf: 0)", 0.0)
h_s = ask("Wasserstand ueber Pumpeneintritt h_s [m]", 0.3)

print("-- Pumpe (Datenblattpunkte, gleiche Anzahl)")
Qp = ask_list("Q [m^3/h]", [0, 2, 4, 6, 8, 10, 12])
Hp = ask_list("H [m]", [12.5, 12.1, 11.2, 9.9, 8.1, 5.9, 3.4])
etap = ask_list("eta [%]", [0, 22, 38, 47, 50, 45, 35]) / 100.0
NPSHp = ask_list("NPSH_R [m]", [0.6, 0.6, 0.7, 0.8, 1.0, 1.3, 1.7])

if not (Qp.size == Hp.size == etap.size == NPSHp.size):
    sys.exit("Fehler: Pumpendaten muessen gleich viele Werte haben.")
if np.any(np.diff(Qp) <= 0):
    sys.exit("Fehler: Q-Werte der Pumpe muessen steigend sein.")
if Q_soll > Qp.max():
    print(f"WARNUNG: Q_soll ({Q_soll}) liegt ausserhalb der Pumpendaten (max {Qp.max()}).")

# ---------- Konstanten ----------
rho, nu, pv, patm, g = 998.0, 1.004e-6, 2339.0, 101325.0, 9.81


def hloss(Qm3h, zE, zV):
    """Verlusthoehe [m] nach Darcy-Weisbach (Swamee-Jain, laminar < 2300)."""
    Qs = np.asarray(Qm3h, dtype=float) / 3600.0
    v = Qs / (np.pi * D**2 / 4.0)
    Re = np.maximum(v * D / nu, 1e-6)
    f_t = 0.25 / np.log10(eps_r / (3.7 * D) + 5.74 / Re**0.9) ** 2
    f = np.where(Re < 2300.0, 64.0 / Re, f_t)
    return (f * L / D + zE + zV) * v**2 / (2.0 * g)


# ---------- Berechnung ----------
Q = np.linspace(0.1, Qp.max(), 400)                  # m^3/h
HP = np.interp(Q, Qp, Hp)
etaP = np.interp(Q, Qp, etap)
NPSHr = np.interp(Q, Qp, NPSHp)
HL = hloss(Q, zeta_E, zeta_V)
HA = H_geo + HL                                      # Anlagenkennlinie
NPSHa = (patm - pv) / (rho * g) + h_s - 0.5 * HL     # halbe Verluste saugseitig

# Betriebspunkt (Schnittpunkt HP = HA)
d = HP - HA
idx = np.where(d[:-1] * d[1:] <= 0)[0]
if idx.size == 0:
    Q_op = H_op = eta_op = NPSHr_op = NPSHa_op = np.nan
else:
    i = idx[0]
    s = d[i] / (d[i] - d[i + 1])
    Q_op = Q[i] + s * (Q[i + 1] - Q[i])
    H_op = np.interp(Q_op, Q, HP)
    eta_op = np.interp(Q_op, Q, etaP)
    NPSHr_op = np.interp(Q_op, Q, NPSHr)
    NPSHa_op = np.interp(Q_op, Q, NPSHa)
P_el = rho * g * (Q_op / 3600.0) * H_op / eta_op     # W

# Drosselstellung, die Q_soll einstellt
A = np.pi * D**2 / 4.0
vS = (Q_soll / 3600.0) / A
HPs = np.interp(Q_soll, Qp, Hp)
Hf0 = hloss(Q_soll, zeta_E, 0.0)
zeta_req = (HPs - H_geo - Hf0) / (vS**2 / (2.0 * g))
HA_req = H_geo + hloss(Q, zeta_E, max(zeta_req, 0.0))

# ---------- Urteil ----------
ok1 = zeta_req >= 0.0                 # Pumpe liefert genug Foerderhoehe
ok2 = zeta_req <= zeta_Vmax           # Pumpe nicht stark ueberdimensioniert
ok3 = NPSHa_op >= 1.3 * NPSHr_op      # NPSH-Reserve
ok4 = eta_op >= 0.7 * etap.max()      # Wirkungsgrad im guten Bereich
verdict = "PUMPE PASST" if all([ok1, ok2, ok3, ok4]) else "PUMPE PASST NICHT"

# ---------- Ausgabe ----------
v_soll = vS
print("\n=== ERGEBNIS ===")
print(f"Geschwindigkeit bei Q_soll:  v = {v_soll:.2f} m/s"
      + ("  (WARNUNG: > 2 m/s)" if v_soll > 2.0 else ""))
print(f"Betriebspunkt (Ventil aktuell): Q = {Q_op:.2f} m^3/h, H = {H_op:.2f} m")
print(f"Pumpenhoehe bei Q_soll:      H = {HPs:.2f} m")
print(f"Noetige Drosselung:          zeta_V = {zeta_req:.2f} (max. {zeta_Vmax:.1f})")
print(f"NPSH_A = {NPSHa_op:.2f} m, NPSH_R = {NPSHr_op:.2f} m, "
      f"Reserve = {NPSHa_op / NPSHr_op:.2f} (Ziel >= 1.3)")
print(f"Wirkungsgrad = {100 * eta_op:.1f} %,  P_el = {P_el:.0f} W")
print("--- Pruefung ---")
print(f"Foerderhoehe ausreichend:   {'JA' if ok1 else 'NEIN'}")
print(f"nicht ueberdimensioniert:   {'JA' if ok2 else 'NEIN'}")
print(f"NPSH-Reserve:               {'JA' if ok3 else 'NEIN'}")
print(f"Wirkungsgrad gut:           {'JA' if ok4 else 'NEIN'}")
print(f"Urteil: {verdict}")

# ---------- Plot ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.plot(Q, HP, "b-", lw=2, label="Pumpe H_P(Q)")
ax1.plot(Q, HA, "r--", lw=1.5, label="Anlage (Ventil aktuell)")
ax1.plot(Q, HA_req, "m:", lw=1.8, label="Anlage (Drossel auf Q_soll)")
ax1.plot(Qp, Hp, "bo", label="Datenblatt")
ax1.plot(Q_op, H_op, "kp", ms=14, mfc="y", label="Betriebspunkt")
ax1.axvline(Q_soll, color="g", lw=1.5, label="Q_soll")
ax1.set(xlabel="Q [m³/h]", ylabel="H [m]", title=f"Urteil: {verdict}")
ax1.grid(True)
ax1.legend(loc="upper right", fontsize=8)

ax2.plot(Q, 100 * etaP, "g-", lw=1.8, label="η")
ax2.set(xlabel="Q [m³/h]", ylabel="η [%]", title="Wirkungsgrad & NPSH")
ax2.grid(True)
ax2b = ax2.twinx()
ax2b.plot(Q, NPSHr, "k-", lw=1.5, label="NPSH_R")
ax2b.plot(Q, NPSHa, "k--", lw=1.5, label="NPSH_A")
ax2b.set_ylabel("NPSH [m]")
h1, l1 = ax2.get_legend_handles_labels()
h2, l2 = ax2b.get_legend_handles_labels()
ax2.legend(h1 + h2, l1 + l2, loc="best", fontsize=8)

plt.tight_layout()
plt.savefig("pumpenauswahl.png", dpi=150)
print("\nPlot gespeichert: pumpenauswahl.png")
plt.show()