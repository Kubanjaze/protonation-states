# Phase 32 — Protonation State Enumerator

**Version:** 1.1 | **Tier:** Micro | **Date:** 2026-03-26

## Goal
Enumerate dominant protonation states for each compound at physiological pH (7.4).
Report pKa estimates, dominant form at pH 7.4, charge state, and visualize
the charge distribution across the library.

CLI: `python main.py --input data/compounds.csv --ph 7.4`

Outputs: protonation_states.csv, charge_distribution.png

## Logic
- Use RDKit's `MolStandardize` (`Normalizer`, `LargestFragmentChooser`, `Uncharger`)
- Identify ionizable groups via SMARTS: amines, carboxylic acids, phenols, sulfonamides, imidazoles
- Estimate pKa using group-specific defaults (literature values):
  - Aliphatic amine: pKa = 10.0 (protonated at pH 7.4 → positive)
  - Aromatic amine: pKa = 4.5 (neutral at pH 7.4)
  - Carboxylic acid: pKa = 4.5 (deprotonated at pH 7.4 → negative)
  - Phenol: pKa = 10.0 (neutral at pH 7.4)
  - Sulfonamide NH: pKa = 10.5 (neutral at pH 7.4)
  - Imidazole: pKa = 6.0 (partially protonated at pH 7.4)
- Formal charge at pH 7.4 = sum of +1 for protonated amines, -1 for deprotonated acids
- Majority rule: if pKa > pH → protonated; if pKa < pH → deprotonated

## Outputs
- protonation_states.csv: compound_name, family, n_ionizable_groups, formal_charge_ph74, dominant_form
- charge_distribution.png: stacked bar or histogram of charge states (−2, −1, 0, +1, +2) by family

## Key Concepts
- RDKit MolStandardize (Normalizer, LargestFragmentChooser, Uncharger)
- SMARTS-based ionizable group detection (amines, acids, phenols, sulfonamides, imidazoles)
- Group-specific pKa estimation using literature default values
- Henderson-Hasselbalch logic: pKa > pH = protonated, pKa < pH = deprotonated

## Verification Checklist
- [x] 45/45 compounds parsed without SMILES failures
- [x] protonation_states.csv contains per-compound ionizable group counts and formal charges
- [x] charge_distribution.png saved to output/
- [x] All compounds are neutral at pH 7.4 (formal charge = 0)
- [x] bzim family correctly identified as having 1 imidazole NH per compound

## Risks
- pKa values are group-specific defaults, not compound-specific predictions (real pKa shifts with environment)
- Only 6 ionizable group types covered; rare functional groups may be missed

## Actual Results (v1.1)

| Family | Mean ionizable groups | Min | Max |
|---|---|---|---|
| benz | 0.08 | 0 | 1 |
| bzim | 1.00 | 1 | 1 |
| ind | 0.00 | 0 | 0 |
| naph | 0.00 | 0 | 0 |
| pyr | 0.00 | 0 | 0 |
| quin | 0.00 | 0 | 0 |

**All 45 compounds: formal charge 0 at pH 7.4**
**Key insight:** This CETP library is predominantly neutral at physiological pH. bzim family has 1 imidazole NH per compound (pKa=6.0 < 7.4 → neutral form dominates). One benz compound has an amine substituent.
**45/45 valid** — no SMILES failures
