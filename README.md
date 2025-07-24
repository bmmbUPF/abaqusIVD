# abaqusIVD

Coupled **mechano‑transport** finite‑element subroutines for the intervertebral disc (IVD).

The repository contains:

- **`Sub_MechDisc.f`** — UMAT for the poro–hyper–visco–elastic mechanical baseline  
- **`Sub_TransDisc.f`** — UEL for fully coupled, multi‑species (O₂, glucose, lactate) transport  
- **`GENERIC_Mechanic.inp` / `GENERIC_Transport.inp`** — Generic morphology example to run the pipeline end‑to‑end

---

## Please cite

If you use any part of these subroutines or scripts, **please cite the paper describing the method**:

1. **Mechanical model**  
   Muñoz‑Moya E., Rasouligandomani M., Ruiz Wills C., Chemorion FK., Piella G., Noailly J. (2024).  
   *Unveiling interactions between intervertebral disc morphologies and mechanical behavior through personalized finite element modeling.*  
   **Front. Bioeng. Biotechnol. 12:1384599.** https://doi.org/10.3389/fbioe.2024.1384599

2. **Transport model (under review)**  
   Muñoz‑Moya, E., Ruiz Wills, C., Nguyen, H. H., Tiulpin, A., Piella, G., & Noailly, J. (2025).  
   *Novel finite element fully coupled multi‑species mechano‑transport simulations approach of the Intervertebral Disc.*  
   **Computer Methods and Programs in Biomedicine** (under review).

If you **use any of the morphed, patient‑personalized FE models from SpineView**, please also cite:

1. **Muñoz‑Moya E, Rasouligandomani M, Ruiz Wills C, Chemorion FK, Piella G, Noailly J (2024).**  
   *Unveiling interactions between intervertebral disc morphologies and mechanical behavior through personalized finite element modeling.*  
   **Front. Bioeng. Biotechnol. 12:1384599.** https://doi.org/10.3389/fbioe.2024.1384599

2. **Muñoz‑Moya E., Rasouligandomani M., Ruiz Wills C., Chemorion F., Piella G., Noailly J. (2023).**  
   *Repository of IVD Patient‑Specific FE Models.* Zenodo. https://doi.org/10.5281/zenodo.8325042

Explore the cohort online at **SpineView UI**: <https://ivd.spineview.upf.edu/>

---

## Repository layout

```
abaqusIVD/
├── GENERIC_Mechanic.inp      # Mechanical baseline (UMAT)
├── GENERIC_Transport.inp     # Multi-species transport (UEL, submodeling)
├── Sub_MechDisc.f            # UMAT
├── Sub_TransDisc.f           # UEL
├── src/                      # (optional) helper scripts/macros
└── README.md
```

---

## Prerequisites

- **Abaqus/Standard 2020 (or newer)**.
- A **Fortran compiler** compatible with your Abaqus version (typically GFortran, or Intel ifort).
- A machine with multiple CPUs if you want to leverage `cpus=Ncpu` and `mp_mode=THREADS`.

---

## Quick start (GENERIC example)

### 1) Run the mechanical baseline (UMAT)

From the repository root:

```bash
abaqus job=GENERIC_Mechanic input=GENERIC_Mechanic user=Sub_MechDisc.f cpus=Ncpu mp_mode=THREADS interactive
```

This creates `GENERIC_Mechanic.odb`, which is then used as the **global model** for transport.

### 2) Run the coupled multi‑species transport (UEL + submodeling)

```bash
abaqus job=GENERIC_Transport input=GENERIC_Transport globalmode=GENERIC_Mechanic.odb user=Sub_TransDisc.f cpus=Ncpu mp_mode=THREADS interactive
```

---

## Reproducing the paper’s workflow (very short version)

1. Run the **mechanical** simulation once to get the desired mechanical baseline.  
2. Launch the **transport** simulation(s), pointing to that mechanical `.odb` with `globalmode=` to reuse the deformations while freely exploring transport parameters and time stepping.  
3. (Optional) Map the solute fields from the GENERIC model to any **SpineView** patient‑personalized model (shared mesh topology enables a direct mapping).

---

## Contact

Questions, issues, or suggestions: **estefano.munoz.moya@gmail.com** or **estefano.munoz@upf.edu**

---

## Funding & acknowledgements

This work received support from:
- EU Horizon 2020 Marie Skłodowska‑Curie (MSCA‑2020‑ITN‑ETN GA: 955735)  
- European Research Council (ERC‑2021‑CoG‑O‑Health‑101044828): Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

---

## Changelog

- **v0.1** – Initial public release with GENERIC example, UMAT (`Sub_MechDisc.f`) and UEL (`Sub_TransDisc.f`).
