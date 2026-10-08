# Pumpen-Simulation

> **ACHTUNG: Alle Werte in `Pumpe_Daten.csv` sind FAKE-Testdaten (Platzhalter/Annahmen), nur zum Testen, ob die Simulation läuft. Vor echter Auswertung durch Datenblattwerte ersetzen.**

## Ordnerstruktur

```
Simulation/
├── README.md
├── Pumpe_Daten.csv        # gemeinsame Eingabe für MATLAB und NumPy
├── Matlab/
│   ├── pumpe_sim.m        # MATLAB-Simulation
│   ├── Ergebnisse.csv     # wird beim Lauf erzeugt
│   └── Plots/             # wird beim Lauf erzeugt
└── NumPy/
    ├── pumpe_sim.py       # NumPy-Simulation
    ├── Ergebnisse.csv     # wird beim Lauf erzeugt
    └── Plots/             # wird beim Lauf erzeugt
```

## Eingabe: `Pumpe_Daten.csv`

Spalten: `Parameter, Wert, Einheit, Beschreibung, Quelle`

Nur Zeilen mit Zahlen in `Wert` werden verwendet. Textzeilen (Typ, Werkstoff, Medium) sind nur Info.

Werte ändern, Datei speichern, Simulation neu starten.

## Ausführen

**MATLAB** (R2020a oder neuer):
1. Ordner `Matlab` öffnen
2. `pumpe_sim` ausführen

**NumPy** (benötigt `numpy` und `matplotlib`):
```
cd NumPy
python pumpe_sim.py
```

## Ausgabe

- `Ergebnisse.csv`: Arbeitspunkt, Leistung, Wirkungsgrad, Anlaufzeit, Temperatur
- `Plots/`: 5 PNG-Dateien (Kennlinien, Wirkungsgrad/Leistung, Ventilstellung, Anlauf, Temperatur)

## Modell in Kürze

- Pumpenkennlinie H(Q) = H0·(1−(Q/Q0)²), Anlagenkennlinie H = H_geo + k(θ)·Q²
- Arbeitspunkt geschlossen berechnet, θ = Drosselventilstellung
- Anlauf: explizites Euler, dt = 0,1 ms
- Temperatur: 1 h, Hysterese-Heizen/Kühlen, Pumpenverlust geht ins Wasser

## Hinweise

- MATLAB-Version ist ungetestet (kein Matlab in der Entwicklungsumgebung). NumPy-Version ist getestet.