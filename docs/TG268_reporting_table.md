# Monte Carlo reporting table (AAPM TG-268 / RECORDS)

Fill this in as values become final. It becomes Table 1 of Paper 1.
Reference: Sechopoulos et al., *Med Phys* 2018;45:e1–e5, <https://doi.org/10.1002/mp.12702>

| Item | Description | Value used | References |
|---|---|---|---|
| Code, version | Monte Carlo code and version | OpenTOPAS v4.3.1 / Geant4 11.4.3 (confirm) | Perl 2012; Faddegon 2020 |
| Validation | How the code/model was validated | Water-tank PDDs, profiles, OFs vs reference (this work) | |
| Timing | CPU time per history / per plan | | |
| Geometry | Linac components modelled; dimensions and materials | Jaws (TsJaws), MLC (TsDivergingMLC, M120), water tank 50×50×40 cm³ | Varian MC package |
| Source | Phase space: origin, plane, number of particles, recycling | Varian TrueBeam 6X IAEA phase space, plane at __ cm, no recycling (chunks) | |
| Cross-sections | Physics data / models | g4em-standard_opt4 | |
| Transport parameters | Production cuts, step limits | 0.1 mm range cut (test 0.05 mm) | |
| Variance reduction | Techniques | none in water validation | |
| Scored quantities | Quantity, scoring geometry, voxel size | DoseToMedium in water; PDD column, profile strips, OF voxel (see config) | |
| Number of histories | Per field / per plan | | |
| Statistical uncertainty | Method and achieved values | Batch method over independent jobs; median __ % | |
| Post-processing | Normalisation, smoothing, fitting | PDD normalised at dmax (local cubic fit); profiles to CAX | |
