# Figure 1A Challenge Diagram

## Script
- `scripts/draw_figure1a_challenge_svg.py`

## Outputs
- `results/Figure1A/Figure1A_challenge_diagram.svg`
- `results/Figure1A/Figure1A_challenge_diagram.pdf`
- `results/Figure1A/Figure1A_challenge_diagram.png`
- `results/Figure1A/README.md`

## Re-run
From the project root:

```bash
python scripts/draw_figure1a_challenge_svg.py
```

The script writes a native SVG first with `svgwrite` when available, then uses `rsvg-convert` to export PDF and PNG. The SVG is composed of independent text, rounded-rectangle, node, edge, legend, and arrow elements so it can be edited in vector software.

## How To Edit
- Colors: edit the `COLORS` dictionary near the top of `scripts/draw_figure1a_challenge_svg.py`.
- Main titles and panel titles: edit `draw_title()`, `draw_left_panel()`, `draw_gap()`, and `draw_right_panel()`.
- Node layout: edit the `local_nodes` and `edges` lists in `draw_left_panel()` and `draw_right_panel()`.
- Canvas size and spacing: edit `WIDTH`, `HEIGHT`, and the panel coordinates inside the drawing functions.

## Design Notes
This is the earlier Figure 1A challenge concept panel: a model yeast essential-gene network, a transfer / annotation gap, and a Fusarium genome-wide prediction network.
