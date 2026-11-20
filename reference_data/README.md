# Reference beam data

Put measured or published reference data here, one folder per beam (e.g. `6X/`, `6FFF/`). These files are **git-ignored by default**: Golden Beam Data and clinic measurements usually cannot be redistributed. Commit only data you are allowed to share.

All files are CSV with a header row. Doses can be in any unit, because the analysis re-normalises them.

| File | Columns | Example name |
|---|---|---|
| PDD | `depth_cm,dose` | `pdd_jaw10cm.csv` |
| Crossline profile | `position_cm,dose` | `profile_x_jaw10cm_d10cm.csv` |
| Inline profile | `position_cm,dose` | `profile_y_jaw10cm_d1p5cm.csv` |
| Output factors | `field,output_factor` (relative to `jaw10cm`) | `output_factors.csv` |

- **Field names** match the validation field names: `jaw1cm`, `jaw2cm`, ..., `jaw30cm`, `mlc10cm`, and so on.
- **Depths** are written as `1p5cm`, `5cm`, `10cm`, `20cm`.
- **Conditions:** the reference must be measured at SSD 100 cm, with depths and output-factor depth (10 cm) matching `config/validation_water.yaml`. If your data uses other conditions, change the config to match.
- **Small fields:** output factors for fields ≤ 3 cm must be detector-corrected *field output factors* (IAEA TRS-483), not raw detector readings.
