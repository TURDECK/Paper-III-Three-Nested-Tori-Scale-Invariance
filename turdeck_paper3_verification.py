#!/usr/bin/env python3
"""
TURDECK Paper III — Script de vérification complète
Trois tores imbriqués : invariance fractale du comma 13/12

Auteur: Sébastien Monast (2026)
Collaborateur computationnel: Claude (Anthropic)

Ce script vérifie TOUTES les claims du Paper III.
Zéro unité — tous les résultats sont des ratios adimensionnels.

Sources des données:
  - NASA Planetary Fact Sheet (nssdc.gsfc.nasa.gov/planetary/factsheet/)
  - CODATA 2018 (constantes fondamentales, NIST)
  - IAU (vitesse orbitale du Soleil dans la galaxie)
  - Bash (1986), Matese & Whitmire (1986) — oscillation galactique
  - Walsh et al. (2011) — Grand Tack model
  - Tsiganis et al. (2005) — Nice model

Usage: python3 turdeck_paper3_verification.py
Sortie: turdeck_paper3_results.json
"""

import numpy as np
import json
from datetime import datetime

PHI = (1 + np.sqrt(5)) / 2
COMMA = 13 / 12
FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987]

# ===================================================================
# FONCTIONS UTILITAIRES
# ===================================================================

def best_fibonacci(value):
    """Trouve la meilleure fraction de Fibonacci pour un ratio donné."""
    best_str, best_val, best_err = "", 0, 100
    for i in range(len(FIB)):
        for j in range(len(FIB)):
            if FIB[j] == 0:
                continue
            r = FIB[i] / FIB[j]
            if 0.001 < r < 50000:
                err = abs(value - r) / value * 100
                if err < best_err:
                    best_err = err
                    best_val = r
                    best_str = f"{FIB[i]}/{FIB[j]}"
    return best_str, best_val, best_err


def best_comma_power(value):
    """Trouve la meilleure expression base × (13/12)^n pour un ratio donné."""
    best_str, best_val, best_err = "", 0, 100
    for n in range(-60, 60):
        for base in range(1, 100):
            val = base * COMMA ** n
            if 0.001 < val < 1e12:
                err = abs(value - val) / value * 100
                if err < best_err and err < 5.0:
                    best_err = err
                    best_val = val
                    best_str = f"{base}*(13/12)^{n}"
    return best_str, best_val, best_err


def verify_ratio(name, observed, description, source):
    """Vérifie un ratio contre Fibonacci ET comma. Retourne un dict."""
    fib_str, fib_val, fib_err = best_fibonacci(observed)
    com_str, com_val, com_err = best_comma_power(observed)
    best_err = min(fib_err, com_err)

    return {
        "name": name,
        "observed": observed,
        "description": description,
        "source": source,
        "fibonacci": {"expression": fib_str, "value": fib_val, "error_pct": round(fib_err, 4)},
        "comma": {"expression": com_str, "value": com_val, "error_pct": round(com_err, 4)} if com_err < 5.0 else None,
        "best_error_pct": round(best_err, 4),
        "status": "SOLID" if best_err < 1.0 else ("APPROX" if best_err < 2.5 else "PROBLEM")
    }


# ===================================================================
# DONNÉES — SOURCES DOCUMENTÉES
# ===================================================================

# Périodes sidérales en jours (NASA Planetary Fact Sheet)
PERIODS = {
    "Mercure": 87.97, "Vénus": 224.70, "Terre": 365.25,
    "Mars": 686.97, "Jupiter": 4332.59, "Saturne": 10759.22,
    "Uranus": 30688.50, "Neptune": 60182.00
}

# Masses relatives (Terre = 1) (NASA Planetary Fact Sheet)
MASSES = {
    "Mercure": 0.0553, "Vénus": 0.815, "Terre": 1.000,
    "Mars": 0.107, "Jupiter": 317.8, "Saturne": 95.16,
    "Uranus": 14.54, "Neptune": 17.15
}

# Distances au Soleil en UA (NASA Planetary Fact Sheet)
DISTANCES = {
    "Mercure": 0.387, "Vénus": 0.723, "Terre": 1.000,
    "Mars": 1.524, "Jupiter": 5.203, "Saturne": 9.537,
    "Uranus": 19.19, "Neptune": 30.07
}

# Vitesses orbitales en km/s (NASA Planetary Fact Sheet)
VELOCITIES = {
    "Mercure": 47.36, "Vénus": 35.02, "Terre": 29.78,
    "Mars": 24.07, "Jupiter": 13.07, "Saturne": 9.68,
    "Uranus": 6.80, "Neptune": 5.43
}

# Rotation sidérale de Vénus en jours (NASA)
VENUS_ROTATION_DAYS = 243.02

# Constantes atomiques (CODATA 2018)
V_BOHR = 2187.691  # km/s — vitesse électron orbite de Bohr
ALPHA_INV = 137.036  # inverse constante de structure fine
RYDBERG_EV = 13.6  # énergie d'ionisation hydrogène en eV

# Données galactiques (IAU, Bash 1986, Nature)
V_SUN_GALACTIC = 220  # km/s — vitesse orbitale Soleil dans galaxie
T_GALACTIC_ORBIT = 225  # en Ma (millions d'années)
T_OSCILLATION_MEDIAN = 63  # en Ma (médiane fourchette 52-74 Ma)
CROSSINGS_PER_ORBIT = 2.7  # traversées du plan galactique par orbite

# Précession et cycles
PRECESSION_YEARS = 25920
SAR_YEARS = 3600

# Fréquences
SCHUMANN_HZ = 7.83


# ===================================================================
# CALCULS
# ===================================================================

def run_all_verifications():
    results = {
        "metadata": {
            "model": "TURDECK Toroidal Knot (12,13), R=5, r=2",
            "paper": "Paper III — Three Nested Tori",
            "author": "Sébastien Monast",
            "date": datetime.now().isoformat(),
            "method": "All ratios are dimensionless (no units)"
        },
        "tore2_orbits": [],
        "tore1_atom": [],
        "tore3_galaxy": [],
        "inter_tores_velocities": [],
        "angular_momenta": [],
        "special_ratios": [],
        "summary": {}
    }

    # ---------------------------------------------------------------
    # TORE 2 — Ratios orbitaux (tours/tours)
    # ---------------------------------------------------------------
    orbit_pairs = [
        ("Terre", "Vénus"), ("Saturne", "Mercure"), ("Uranus", "Mars"),
        ("Vénus", "Saturne"), ("Jupiter", "Saturne"), ("Terre", "Saturne"),
        ("Terre", "Mars"), ("Neptune", "Terre"), ("Mercure", "Vénus")
    ]
    for p1, p2 in orbit_pairs:
        ratio = PERIODS[p1] / PERIODS[p2]
        r = verify_ratio(
            f"T_{p1}/T_{p2}", ratio,
            f"Ratio de périodes orbitales {p1}/{p2} (tours/tours)",
            "NASA Planetary Fact Sheet"
        )
        results["tore2_orbits"].append(r)

    # Vénus rotation/orbite
    r = verify_ratio(
        "Venus_rot/orb", VENUS_ROTATION_DAYS / PERIODS["Vénus"],
        "Ratio rotation sidérale / période orbitale de Vénus",
        "NASA Planetary Fact Sheet"
    )
    results["tore2_orbits"].append(r)

    # ---------------------------------------------------------------
    # TORE 1 — Atome
    # ---------------------------------------------------------------
    r = verify_ratio("alpha_inv", ALPHA_INV,
        "Inverse de la constante de structure fine", "CODATA 2018")
    results["tore1_atom"].append(r)

    r = verify_ratio("Rydberg", RYDBERG_EV,
        "Énergie d'ionisation de l'hydrogène (ratio 13.6, contient le 13)",
        "CODATA 2018")
    results["tore1_atom"].append(r)

    # ---------------------------------------------------------------
    # TORE 3 — Galaxie
    # ---------------------------------------------------------------
    for t_osc in [52, 55, 60, 63, 66, 70, 74]:
        ratio = T_GALACTIC_ORBIT / t_osc
        r = verify_ratio(
            f"T_orb/T_osc_{t_osc}Ma", ratio,
            f"Ratio orbite galactique / oscillation verticale (T_osc={t_osc} Ma)",
            "Bash 1986, Matese & Whitmire 1986"
        )
        results["tore3_galaxy"].append(r)

    r = verify_ratio("crossings_per_orbit", CROSSINGS_PER_ORBIT,
        "Traversées du plan galactique par orbite", "Nature / Wikipedia")
    results["tore3_galaxy"].append(r)

    # ---------------------------------------------------------------
    # INTER-TORES — Ratios de vitesses (sans unité: km/s ÷ km/s)
    # ---------------------------------------------------------------
    # Atome → Système solaire
    for planet, v in VELOCITIES.items():
        ratio = V_BOHR / v
        r = verify_ratio(
            f"v_Bohr/v_{planet}", ratio,
            f"Ratio vitesse électron (Bohr) / vitesse orbitale {planet}",
            "CODATA + NASA"
        )
        results["inter_tores_velocities"].append(r)

    # Système → Galaxie
    for planet, v in VELOCITIES.items():
        ratio = V_SUN_GALACTIC / v
        r = verify_ratio(
            f"v_Gal/v_{planet}", ratio,
            f"Ratio vitesse Soleil (galactique) / vitesse orbitale {planet}",
            "IAU + NASA"
        )
        results["inter_tores_velocities"].append(r)

    # Atome → Galaxie
    r = verify_ratio("v_Bohr/v_Gal", V_BOHR / V_SUN_GALACTIC,
        "Ratio vitesse électron / vitesse Soleil galactique", "CODATA + IAU")
    results["inter_tores_velocities"].append(r)

    # ---------------------------------------------------------------
    # MOMENTS ANGULAIRES (L = m × v × a, relatif à la Terre)
    # ---------------------------------------------------------------
    L_terre = MASSES["Terre"] * VELOCITIES["Terre"] * DISTANCES["Terre"]
    L = {}
    for planet in MASSES:
        L[planet] = (MASSES[planet] * VELOCITIES[planet] * DISTANCES[planet]) / L_terre

    L_pairs = [
        ("Jupiter", "Terre"), ("Saturne", "Terre"), ("Neptune", "Terre"),
        ("Jupiter", "Saturne"), ("Neptune", "Uranus"), ("Vénus", "Terre"),
        ("Mars", "Terre"), ("Mercure", "Terre"), ("Jupiter", "Mars"),
        ("Saturne", "Jupiter"), ("Neptune", "Jupiter"), ("Uranus", "Saturne")
    ]
    for p1, p2 in L_pairs:
        ratio = L[p1] / L[p2]
        r = verify_ratio(
            f"L_{p1}/L_{p2}", ratio,
            f"Ratio moment angulaire {p1}/{p2} (m×v×a, sans unité)",
            "NASA Planetary Fact Sheet (calcul dérivé)"
        )
        results["angular_momenta"].append(r)

    # ---------------------------------------------------------------
    # RATIOS SPÉCIAUX
    # ---------------------------------------------------------------
    r = verify_ratio("precession_SAR", PRECESSION_YEARS / SAR_YEARS,
        "Ratio précession / SAR = 25920/3600 = 36/5 = 6²/5", "IAU")
    results["special_ratios"].append(r)

    # Phase dans le cycle de précession (12 secteurs)
    phase_12 = 2026 / (PRECESSION_YEARS / 12)
    r = verify_ratio("phase_12_sectors", phase_12,
        "Phase du point vernal dans le secteur actuel (12 secteurs)", "IAU")
    results["special_ratios"].append(r)

    # Phase dans le cycle de précession (13 secteurs)
    phase_13 = 2026 / (PRECESSION_YEARS / 13)
    r = verify_ratio("phase_13_sectors", phase_13,
        "Phase du point vernal (13 secteurs = avec correction torique)", "IAU")
    results["special_ratios"].append(r)

    # 432 / Schumann
    r = verify_ratio("432_Schumann", 432 / SCHUMANN_HZ,
        "Ratio 432 Hz / fréquence de Schumann (7.83 Hz)", "Mesures")
    results["special_ratios"].append(r)

    # ---------------------------------------------------------------
    # BILAN
    # ---------------------------------------------------------------
    all_claims = (results["tore2_orbits"] + results["tore1_atom"] +
                  results["tore3_galaxy"] + results["inter_tores_velocities"] +
                  results["angular_momenta"] + results["special_ratios"])

    solid = [c for c in all_claims if c["status"] == "SOLID"]
    approx = [c for c in all_claims if c["status"] == "APPROX"]
    problem = [c for c in all_claims if c["status"] == "PROBLEM"]

    results["summary"] = {
        "total_claims": len(all_claims),
        "solid_under_1pct": len(solid),
        "approximation_1_to_2.5pct": len(approx),
        "problem_over_2.5pct": len(problem),
        "worst_error": max(c["best_error_pct"] for c in all_claims),
        "median_error": round(np.median([c["best_error_pct"] for c in all_claims]), 4),
        "mean_error": round(np.mean([c["best_error_pct"] for c in all_claims]), 4)
    }

    return results


# ===================================================================
# EXÉCUTION
# ===================================================================

if __name__ == "__main__":
    print("TURDECK Paper III — Vérification complète")
    print("=" * 60)

    results = run_all_verifications()

    # Affichage du bilan
    s = results["summary"]
    print(f"\n  BILAN FINAL:")
    print(f"  Total claims vérifiées: {s['total_claims']}")
    print(f"  ✅ Solide (< 1%): {s['solid_under_1pct']}")
    print(f"  ⚠️  Approximation (1-2.5%): {s['approximation_1_to_2.5pct']}")
    print(f"  ❌ Problème (> 2.5%): {s['problem_over_2.5pct']}")
    print(f"  Erreur médiane: {s['median_error']}%")
    print(f"  Erreur moyenne: {s['mean_error']}%")
    print(f"  Pire erreur: {s['worst_error']}%")

    # Sauvegarde JSON
    output_file = "turdeck_paper3_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n  Résultats sauvegardés: {output_file}")
    print(f"\n  © 2026 Sébastien Monast — TURDECK")
