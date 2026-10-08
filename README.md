# Cooling-System 2026

Kühlsystem-Prüfstand mit Sensorik und Messdatenerfassung.

## Aufbau

Becken → Rohrleitung → Pumpe → Druckmessung

## Komponenten

- **Tank / Becken**
- **Pumpe** (korrosionsbeständig)
- **Regelventil / Drosselventil**
- **Heizelement**

## Sensorik

| Sensor | Anzahl / Position |
|---|---|
| Volumenstromsensor | 1 |
| Drucksensor | 1 vor und 1 nach der Pumpe |
| Temperatursensor | Einlass und Auslass |
| Elektrische Messtechnik | Strom und Spannung |

## Messtechnik / Datenerfassung

Auswahl je nach Budget:

- **Arduino / Raspberry Pi** (günstig)
- **DAQ** (Data Acquisition, teuer, präzise)

## Datenverarbeitung

- **MATLAB** (kompatibel mit vielen Geräten und DAQ-Systemen)

## Offene Punkte

- Abgleich mit alter Pumpe bzw. neuem Inverter
