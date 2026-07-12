# Figure 1B Reconstruction Diagram

## Script
- `scripts/draw_figure1b_reconstruction_svg.py`

## Outputs
- `results/Figure1B/Figure1B_reconstruction_diagram.svg`
- `results/Figure1B/Figure1B_reconstruction_diagram.pdf`
- `results/Figure1B/Figure1B_reconstruction_diagram.png`
- `results/Figure1B/README.md`

## Re-run
From the project root:

```bash
python scripts/draw_figure1b_reconstruction_svg.py
```

The script writes a native SVG first with `svgwrite` when available, then uses `rsvg-convert` to export PDF and PNG. The SVG is composed of independent text, rounded-rectangle, node, edge, legend, and arrow elements so it can be edited in vector software. The generated PDF keeps embedded TrueType fonts when checked with `pdffonts`.

## How To Edit
- Colors: edit the `COLORS` dictionary near the top of `scripts/draw_figure1b_reconstruction_svg.py`.
- Main titles and panel title: edit `draw_title()`.
- Node text and layout: edit `draw_main_chain()` and `draw_phi_branch()`.
- Canvas size and spacing: edit `WIDTH`, `HEIGHT`, and the panel coordinates inside the drawing functions.

## Design Notes
The diagram is an editorial-style scientific figure panel showing a single left-to-right reconstruction workflow: yeast essential evidence passes through orthogroup transfer and canonical Fusarium mapping, PHI-base lethal evidence joins at the aggregation step, and the workflow terminates in one reconstructed essential gene set.
