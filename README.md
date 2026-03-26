# protonation-states — Phase 32

Enumerate dominant protonation states at physiological pH (7.4) using rule-based pKa estimates.
Reports ionizable groups, formal charge per compound, and charge distribution by scaffold family.

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv --ph 7.4
```

## Outputs

| File | Description |
|---|---|
| `output/protonation_states.csv` | Compound × ionizable groups + formal charge at pH 7.4 |
| `output/charge_distribution.png` | Stacked bar: charge states by scaffold family |

## pKa Rules

| Group | pKa | Form at pH 7.4 |
|---|---|---|
| Aliphatic amine | 10.0 | Protonated (+1) |
| Aromatic amine | 4.5 | Neutral |
| Carboxylic acid | 4.5 | Deprotonated (−1) |
| Phenol | 10.0 | Neutral |
| Sulfonamide NH | 10.5 | Neutral |
| Imidazole | 6.0 | Neutral (pKa < 7.4) |
