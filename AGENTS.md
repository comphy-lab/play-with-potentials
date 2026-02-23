# play-with-potentials

Python scripts for plotting and comparing effective interaction potentials.

## Structure

- `plot_u.py`: Main two-panel plot (vary d and rho).
- `plot_u_v1.py`: Variant sweep over d.
- `plot_u_v0.py`: Baseline single-case plot.
- `u_of_r.pdf`: Generated figure output.

## Development

```bash
# Environment
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib

# Run main script
python3 plot_u.py
```

## Guidelines

- Style: Follow PEP 8 for Python changes.
- Math plotting: Keep formulas numerically safe near singular points (`r > 2d` cutoffs).
- Output: Script outputs should remain reproducible and write figures to repository root unless changed deliberately.
