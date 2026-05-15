[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Tight low complexity approximation for the $Qm||C_{max}$ problem via mathematical programming modeling

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Tight low complexity approximation for the $Qm||C_{max}$ problem via mathematical programming modeling](https://doi.org/10.1287/ijoc.2025.1129) by L. Savant Aira, R. Scatamacchia, and F. Della Croce.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2025.1129

https://doi.org/10.1287/ijoc.2025.1129.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{TightLowComplexityApproximation,
  author =        {L. {Savant Aira} and R. Scatamacchia and F. Della Croce},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Tight Low Complexity Approximation for the $Qm||C_{max}$ Problem via Mathematical Programming Modeling}},
  year =          {2025},
  doi =           {10.1287/ijoc.2025.1129.cd},
  url =           {https://github.com/INFORMSJoC/2025.1129},
  note =          {Available for download at https://github.com/INFORMSJoC/2025.1129},
}  
```

## Description

This repository contains Gurobi models used to study worst-case approximation ratios for low-complexity algorithms for the uniform machine scheduling problem $Qm||C_{max}$, including LPT, MLPT, FFD, and variants.

## Requirements

- Python >= 3.12
- Gurobi with a valid license. We suggest using Gurobi 13.0.1 with an "Academic WLS License".

Install dependencies with:

```bash
uv sync
```

or:

```bash
pip install gurobipy==13.0.1
```

## Reproduce

To reproduce the results of Table 1 of the paper, run:

```bash
uv run src/LPT.py
```

To reproduce the results of Table 2 of the paper, run:

```bash
uv run src/MLPT.py
```

To reproduce the results of Table 3 of the paper (column by column), run:

```bash
uv run src/FFD.py
uv run src/FFD_reversed.py
uv run src/FFD_bidirectional.py
uv run src/COMBINE.py
uv run src/COMBINE_bidirectional.py
```

Note that some larger cases are commented out in `LPT.py` and `MLPT.py` because they are computationally expensive. As reported in the paper, we used a server-grade Intel Xeon Gold 6130 CPU with 128GB of RAM and 16 cores to solve also these larger cases.

## Structure

- `src/`: mathematical programming models for LPT, MLPT, FFD, and variants
- `AUTHORS`, `LICENSE`: authorship and license information
- `README.md`: this file
