"""
Pumpen-Simulation (NumPy + Matplotlib)
Kennlinien, Arbeitspunkte, Wirkungsgrad, Anlauf (Euler), Temperaturregelung.
Aufruf: python pumpe_sim.py   -> liest ../Pumpe_Daten.csv, schreibt Plots/ und Ergebnisse.csv
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PFAD = os.path.join(HERE, "..", "Pumpe_Daten.csv")
PLOT_DIR = os.path.join(HERE, "Plots")
os.makedirs(PLOT_DIR, exist_ok=True)


# ---------- Parameter ----------
def lade_parameter(pfad):
    p = {}
    with open(pfad, newline="", encoding="utf-8") as f:
        for zeile in csv.DictReader(f):
            try:
                p[zeile["Parameter"]] = float(zeile["Wert"])
            except ValueError:
                pass  # Textwerte (Typ, Werkstoff, Medium) ueberspringen
    return p


P = lade_parameter(CSV_PFAD)
rho, g = P["rho"], P["g"]
H0, Q0, ETA_MAX, Q_BEP = P["H0"], P["Q0"], P["eta_max"], P["Q_bep"]
H_GEO, K_ROHR, K_V0 = P["H_geo"], P["k_rohr"], P["k_ventil0"]
A = np.pi * (P["D_innen"] / 2) ** 2  # Rohrquerschnitt [m2]
L = P["L_rohr"]


# ---------- Modell ----------
def H_pumpe(Q):
    """Pumpenkennlinie H(Q) [m]"""
    return H0 * (1 - (Q / Q0) ** 2)


def eta(Q):
    """Wirkungsgrad [-], Glockenkurve um Bestpunkt"""
    return ETA_MAX * np.exp(-((Q - Q_BEP) / (0.5 * Q_BEP)) ** 2)


def P_hyd(Q, H):
    """Hydraulische Leistung [W]"""
    return rho * g * Q * H


def P_el(Q):
    """Elektrische Aufnahmeleistung [W]"""
    return P_hyd(Q, H_pumpe(Q)) / eta(Q)


def k_gesamt(theta):
    """Anlagenbeiwert [s2/m5] bei Ventilstellung theta (0..1)"""
    return K_ROHR + K_V0 * (1 / theta ** 2 - 1)


def H_anlage(Q, theta):
    """Anlagenkennlinie [m]"""
    return H_GEO + k_gesamt(theta) * Q ** 2


def arbeitspunkt(theta):
    """Schnitt Pumpe/Anlage, geschlossen geloest: Q_op, H_op"""
    Qop = np.sqrt((H0 - H_GEO) / (H0 / Q0 ** 2 + k_gesamt(theta)))
    return Qop, H_pumpe(Qop)


def anlauf(theta=1.0, T_end=0.5, dt=1e-4):
    """Anlauf aus Ruhe: dQ/dt = g*A/L * (H_pumpe - H_anlage), explizites Euler"""
    n = int(round(T_end / dt)) + 1
    t = np.arange(n) * dt
    Q = np.zeros(n)
    faktor = g * A / L
    for i in range(n - 1):
        Q[i + 1] = Q[i] + dt * faktor * (H_pumpe(Q[i]) - H_anlage(Q[i], theta))
    return t, Q


def temperatur(T_end=3600.0, dt=1.0):
    """Beckentemperatur: Hysterese-Regler Heizen/Kuehlen, Pumpenverlust ins Wasser"""
    Qop, Hop = arbeitspunkt(1.0)
    Q_verl = float(P_el(Qop) - P_hyd(Qop, Hop))  # Verlustleistung Pumpe [W]
    C = rho * P["V_tank"] * P["c_w"]              # Waermekapazitaet Becken [J/K]
    ts, hy = P["T_soll"], P["Hysterese"]
    n = int(round(T_end / dt)) + 1
    t = np.arange(n) * dt
    T = np.empty(n)
    T[0] = P["T_start"]
    heiz = np.zeros(n, dtype=bool)
    kuehl = np.zeros(n, dtype=bool)
    h_an = k_an = False
    for i in range(n - 1):
        Ti = T[i]
        if Ti < ts - hy:
            h_an = True
        elif Ti >= ts:
            h_an = False
        if Ti > ts + hy:
            k_an = True
        elif Ti <= ts:
            k_an = False
        q = (P["P_heiz"] if h_an else 0.0) - (P["P_kuehl"] if k_an else 0.0)
        q += Q_verl - P["UA"] * (Ti - P["T_umg"])
        T[i + 1] = Ti + dt * q / C
        heiz[i], kuehl[i] = h_an, k_an
    heiz[-1], kuehl[-1] = h_an, k_an
    return t, T, heiz, kuehl, Q_verl


# ---------- Rechnung ----------
def main():
    Qop, Hop = arbeitspunkt(1.0)
    Phyd_op = float(P_hyd(Qop, Hop))
    Pel_op = float(P_el(Qop))
    eta_op = float(eta(Qop))

    tA, QA = anlauf(1.0)
    t95 = tA[np.argmax(QA >= 0.95 * Qop)]

    tT, TT, heiz, kuehl, Q_verl = temperatur()

    ergebnisse = [
        ("Foerderstrom_Q_op", Qop * 6e4, "L/min"),
        ("Foerderhoehe_H_op", Hop, "m"),
        ("Hydraulische_Leistung", Phyd_op, "W"),
        ("Wirkungsgrad_op", eta_op, "-"),
        ("Elektrische_Leistung_op", Pel_op, "W"),
        ("Verlustleistung_Pumpe_ins_Wasser", Q_verl, "W"),
        ("Anlaufzeit_t95", t95, "s"),
        ("Heizelement_Anteil", heiz.mean() * 100, "%"),
        ("Kuehlung_Anteil", kuehl.mean() * 100, "%"),
        ("Temperatur_Ende", TT[-1], "C"),
        ("Temperatur_Max", TT.max(), "C"),
    ]
    with open(os.path.join(HERE, "Ergebnisse.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Groesse", "Wert", "Einheit"])
        for name, wert, einheit in ergebnisse:
            w.writerow([name, f"{wert:.4f}", einheit])
    for name, wert, einheit in ergebnisse:
        print(f"{name:35s} {wert:12.4f} {einheit}")

    # ---- Plots ----
    Qs = np.linspace(0, Q0, 400)
    thetas = [1.0, 0.75, 0.5, 0.3]

    # 1 Kennlinien + Arbeitspunkte
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Qs * 6e4, H_pumpe(Qs), "k", lw=2, label="Pumpenkennlinie")
    for th in thetas:
        ax.plot(Qs * 6e4, H_anlage(Qs, th), "--", label=f"Anlage theta={th:.2f}")
        Qo, Ho = arbeitspunkt(th)
        ax.plot(Qo * 6e4, Ho, "o", ms=8, label=f"AP theta={th:.2f}: {Qo*6e4:.1f} L/min, {Ho:.1f} m")
    ax.set_xlabel("Foerderstrom Q [L/min]")
    ax.set_ylabel("Foerderhoehe H [m]")
    ax.set_title("Pumpen- und Anlagenkennlinie, Arbeitspunkte")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "1_Kennlinien.png"), dpi=150)
    plt.close(fig)

    # 2 Wirkungsgrad + elektrische Leistung
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    a1.plot(Qs * 6e4, eta(Qs), "g", lw=2)
    a1.axvline(Q_BEP * 6e4, color="gray", ls=":")
    a1.set_ylabel("Wirkungsgrad eta [-]")
    a1.grid(alpha=0.3)
    a2.plot(Qs * 6e4, P_el(Qs), "r", lw=2, label="P_el")
    a2.axhline(P["P_nenn"], color="k", ls="--", label="P_nenn")
    a2.set_xlabel("Foerderstrom Q [L/min]")
    a2.set_ylabel("Elektrische Leistung [W]")
    a2.grid(alpha=0.3)
    a2.legend()
    fig.suptitle("Wirkungsgrad und Leistungsaufnahme")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "2_Wirkungsgrad_Leistung.png"), dpi=150)
    plt.close(fig)

    # 3 Ventilstellung -> Arbeitspunkt
    th = np.linspace(0.2, 1.0, 60)
    Qv, Hv = arbeitspunkt(th)
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    a1.plot(th * 100, Qv * 6e4, "b", lw=2)
    a1.set_ylabel("Q [L/min]")
    a2.plot(th * 100, Hv, "m", lw=2)
    a2.set_ylabel("H [m]")
    a3.plot(th * 100, P_el(Qv), "r", lw=2)
    a3.set_ylabel("P_el [W]")
    a3.set_xlabel("Ventilstellung theta [%]")
    for a in (a1, a2, a3):
        a.grid(alpha=0.3)
    fig.suptitle("Drosselventil: Arbeitspunkt in Abhaengigkeit der Stellung")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "3_Ventilstellung.png"), dpi=150)
    plt.close(fig)

    # 4 Anlauf
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(tA * 1000, QA * 6e4, "b", lw=2, label="Q(t)")
    ax.axhline(Qop * 6e4, color="k", ls="--", label="Q_op")
    ax.axvline(t95 * 1000, color="gray", ls=":", label=f"t95 = {t95*1000:.0f} ms")
    ax.set_xlabel("Zeit t [ms]")
    ax.set_ylabel("Foerderstrom Q [L/min]")
    ax.set_title("Anlauf aus Ruhe (Euler, dt = 0.1 ms)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "4_Anlauf.png"), dpi=150)
    plt.close(fig)

    # 5 Temperatur
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    a1.plot(tT / 60, TT, "r", lw=2, label="T Becken")
    a1.axhline(P["T_soll"], color="k", ls="--", label="T_soll")
    a1.axhspan(P["T_soll"] - P["Hysterese"], P["T_soll"] + P["Hysterese"], color="gray", alpha=0.15)
    a1.set_ylabel("Temperatur [C]")
    a1.grid(alpha=0.3)
    a1.legend()
    a2.fill_between(tT / 60, 0, heiz.astype(float), step="post", color="orange", alpha=0.6, label="Heizen")
    a2.fill_between(tT / 60, 0, -kuehl.astype(float), step="post", color="blue", alpha=0.6, label="Kuehlen")
    a2.set_xlabel("Zeit t [min]")
    a2.set_ylabel("Zustand")
    a2.set_yticks([-1, 0, 1])
    a2.legend()
    a2.grid(alpha=0.3)
    fig.suptitle("Beckentemperatur mit Hysterese-Regler (1 h)")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "5_Temperatur.png"), dpi=150)
    plt.close(fig)

    print("Plots gespeichert in:", PLOT_DIR)


if __name__ == "__main__":
    main()
