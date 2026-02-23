# play-with-potentials

Python scripts for exploring and visualizing effective interaction potentials and comparing them with Lennard-Jones behavior.

## Overview

This repository contains lightweight analysis scripts that compute attractive/repulsive potential terms and generate publication-style plots of \(U(r)\). It is intended for quick exploratory studies and figure generation.

## Features

- Multiple script variants for potential analysis (`plot_u_v0.py`, `plot_u_v1.py`, `plot_u.py`).
- Comparison against Lennard-Jones reference curves.
- PDF figure export (`u_of_r.pdf`) using LaTeX-rendered labels in Matplotlib.

## Requirements

- Python 3.10+
- `numpy`
- `matplotlib`
- A working LaTeX installation (for `matplotlib.rcParams['text.usetex'] = True`)

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib
```

## Usage

```bash
python3 plot_u.py
```

This writes `u_of_r.pdf` to the repository root.

## Repository Contents

- `plot_u.py`: Two-panel visualization varying both \(d\) and \(\rho\).
- `plot_u_v1.py`: Single-panel sweep over \(d\).
- `plot_u_v0.py`: Baseline single-case plot.
- `u_of_r.pdf`: Example generated output.

## License

MIT (recommended). Add a `LICENSE` file if you plan to distribute this repository publicly.
