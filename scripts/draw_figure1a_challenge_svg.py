#!/usr/bin/env python3
"""
脚本名称: draw_figure1a_challenge_svg.py
日期: 2026-05-06
作者: OpenAI/Codex + Jie Huang workflow support
功能描述: 直接生成 Figure 1A challenge concept diagram 原生 SVG，并导出 PDF 和 PNG。
"""

from __future__ import annotations

import html
import shutil
import subprocess
from pathlib import Path
from typing import Protocol


REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "results" / "Figure1A"
SVG_OUT = OUT_DIR / "Figure1A_challenge_diagram.svg"
PDF_OUT = OUT_DIR / "Figure1A_challenge_diagram.pdf"
PNG_OUT = OUT_DIR / "Figure1A_challenge_diagram.png"
README_OUT = OUT_DIR / "README.md"

WIDTH = 1600
HEIGHT = 900

# 统一颜色配置，修改配色优先改这里
COLORS = {
    "ink": "#232323",
    "muted": "#5A5F66",
    "frame": "#C7CDD4",
    "soft_shadow": "#D8DEE7",
    "left_header": "#E9F2FB",
    "right_header": "#FCECEC",
    "essential_blue": "#4F8BCB",
    "pred_red": "#D85F5A",
    "nonessential": "#DADFE5",
    "nonessential_stroke": "#AEB6BF",
    "edge": "#A7B0BA",
    "gap_box": "#FFFFFF",
    "gap_stroke": "#C8CFD8",
    "arrow": "#B8C0CA",
}


class CanvasProtocol(Protocol):
    def add_def(self, text: str) -> None: ...
    def rect(self, x: float, y: float, w: float, h: float, **attrs: object) -> None: ...
    def circle(self, cx: float, cy: float, r: float, **attrs: object) -> None: ...
    def line(self, x1: float, y1: float, x2: float, y2: float, **attrs: object) -> None: ...
    def text(self, text: str, x: float, y: float, **attrs: object) -> None: ...
    def write(self, path: Path) -> None: ...


class SvgCanvas:
    """Small native SVG writer used when svgwrite is unavailable."""

    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        self.elements: list[str] = []
        self.defs: list[str] = []

    def add_def(self, text: str) -> None:
        self.defs.append(text)

    def add(self, text: str) -> None:
        self.elements.append(text)

    def rect(self, x: float, y: float, w: float, h: float, **attrs: object) -> None:
        self.add(f"<rect {format_attrs({'x': x, 'y': y, 'width': w, 'height': h, **attrs})} />")

    def circle(self, cx: float, cy: float, r: float, **attrs: object) -> None:
        self.add(f"<circle {format_attrs({'cx': cx, 'cy': cy, 'r': r, **attrs})} />")

    def line(self, x1: float, y1: float, x2: float, y2: float, **attrs: object) -> None:
        self.add(f"<line {format_attrs({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2, **attrs})} />")

    def text(self, text: str, x: float, y: float, **attrs: object) -> None:
        self.add(f"<text {format_attrs({'x': x, 'y': y, **attrs})}>{html.escape(text)}</text>")

    def write(self, path: Path) -> None:
        defs = f"<defs>{''.join(self.defs)}</defs>" if self.defs else ""
        svg = (
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}" role="img" aria-label="Figure 1A challenge diagram">\n'
            f"{defs}\n"
            + "\n".join(self.elements)
            + "\n</svg>\n"
        )
        path.write_text(svg, encoding="utf-8")


class SvgwriteCanvas:
    """svgwrite-backed SVG writer that keeps text and shapes editable."""

    def __init__(self, width: int, height: int) -> None:
        import svgwrite

        self.dwg = svgwrite.Drawing(
            filename=str(SVG_OUT),
            size=(width, height),
            viewBox=f"0 0 {width} {height}",
            profile="full",
            debug=False,
        )
        self.dwg.attribs["role"] = "img"
        self.dwg.attribs["aria-label"] = "Figure 1A challenge diagram"

    def add_def(self, text: str) -> None:
        marker = self.dwg.marker(
            id="arrowHead",
            insert=(10, 4.5),
            size=(12, 9),
            orient="auto",
            markerUnits="strokeWidth",
        )
        marker.add(self.dwg.path(d="M0,0 L11,4.5 L0,9 z", fill=COLORS["arrow"]))
        self.dwg.defs.add(marker)

    def rect(self, x: float, y: float, w: float, h: float, **attrs: object) -> None:
        self.dwg.add(self.dwg.rect(insert=(x, y), size=(w, h), **normalize_svgwrite_attrs(attrs)))

    def circle(self, cx: float, cy: float, r: float, **attrs: object) -> None:
        self.dwg.add(self.dwg.circle(center=(cx, cy), r=r, **normalize_svgwrite_attrs(attrs)))

    def line(self, x1: float, y1: float, x2: float, y2: float, **attrs: object) -> None:
        self.dwg.add(self.dwg.line(start=(x1, y1), end=(x2, y2), **normalize_svgwrite_attrs(attrs)))

    def text(self, text: str, x: float, y: float, **attrs: object) -> None:
        self.dwg.add(self.dwg.text(text, insert=(x, y), **normalize_svgwrite_attrs(attrs)))

    def write(self, path: Path) -> None:
        self.dwg.filename = str(path)
        self.dwg.save(pretty=True)


def format_attrs(attrs: dict[str, object]) -> str:
    parts = []
    for key, value in attrs.items():
        if value is None:
            continue
        parts.append(f'{key.replace("_", "-")}="{html.escape(str(value), quote=True)}"')
    return " ".join(parts)


def normalize_svgwrite_attrs(attrs: dict[str, object]) -> dict[str, object]:
    return {key.replace("_", "-"): value for key, value in attrs.items() if value is not None}


def make_canvas() -> CanvasProtocol:
    # 优先使用 svgwrite；缺失时使用内置原生 SVG writer
    try:
        import svgwrite  # noqa: F401

        print("[Figure1A] svgwrite is available; using svgwrite backend.")
        return SvgwriteCanvas(WIDTH, HEIGHT)
    except ImportError:
        print("[Figure1A] svgwrite is not installed; using built-in native SVG writer fallback.")
        return SvgCanvas(WIDTH, HEIGHT)


def add_defs(canvas: CanvasProtocol) -> None:
    # 添加统一箭头定义
    canvas.add_def(
        """
        <marker id="arrowHead" markerWidth="12" markerHeight="9" refX="10" refY="4.5"
                orient="auto" markerUnits="strokeWidth">
          <path d="M0,0 L11,4.5 L0,9 z" fill="#B8C0CA"/>
        </marker>
        """
    )


def text_style(size: int, weight: str = "400", fill: str | None = None, anchor: str = "start") -> dict[str, object]:
    # 统一文本样式，保证 SVG/PDF 中文字尽量保持可编辑
    return {
        "font_family": "Arial, Helvetica, sans-serif",
        "font_size": size,
        "font_weight": weight,
        "fill": fill or COLORS["ink"],
        "text_anchor": anchor,
        "dominant_baseline": "middle",
    }


def draw_title(canvas: CanvasProtocol) -> None:
    # 绘制 panel label 和顶部标题
    canvas.text("A", 64, 70, **text_style(42, "700", COLORS["ink"], "start"))
    canvas.text("Challenge of essential gene prediction", WIDTH / 2, 63, **text_style(34, "700", COLORS["ink"], "middle"))
    canvas.text("in a non-model fungal pathogen", WIDTH / 2, 105, **text_style(28, "600", COLORS["muted"], "middle"))


def draw_panel(canvas: CanvasProtocol, x: int, y: int, w: int, h: int, title: str, header_fill: str) -> None:
    # 绘制主框、柔和标题条和边框
    canvas.rect(x + 8, y + 10, w, h, rx=24, ry=24, fill=COLORS["soft_shadow"], opacity=0.35)
    canvas.rect(x, y, w, h, rx=24, ry=24, fill="#FFFFFF", stroke=COLORS["frame"], stroke_width=1.6)
    canvas.rect(x, y, w, 62, rx=24, ry=24, fill=header_fill, stroke="none")
    canvas.rect(x, y + 38, w, 24, fill=header_fill, stroke="none")
    canvas.line(x, y + 62, x + w, y + 62, stroke=COLORS["frame"], stroke_width=1.2)
    canvas.text(title, x + w / 2, y + 32, **text_style(20, "700", COLORS["ink"], "middle"))


def draw_network(
    canvas: CanvasProtocol,
    origin_x: int,
    origin_y: int,
    nodes: list[dict[str, object]],
    edges: list[tuple[int, int]],
) -> None:
    # 绘制手工布局网络连线
    for source, target in edges:
        n1 = nodes[source]
        n2 = nodes[target]
        canvas.line(
            origin_x + n1["x"],
            origin_y + n1["y"],
            origin_x + n2["x"],
            origin_y + n2["y"],
            stroke=COLORS["edge"],
            stroke_width=2.2,
            opacity=0.72,
            stroke_linecap="round",
        )

    # 绘制手工布局网络节点
    for node in nodes:
        canvas.circle(
            origin_x + node["x"],
            origin_y + node["y"],
            node["r"],
            fill=node["fill"],
            stroke=node.get("stroke", "#FFFFFF"),
            stroke_width=2.2,
        )
        if node.get("label"):
            canvas.text(str(node["label"]), origin_x + node["x"], origin_y + node["y"] + 1, **text_style(11, "700", "#FFFFFF", "middle"))


def draw_legend(canvas: CanvasProtocol, x: int, y: int, items: list[tuple[str, str, str]]) -> None:
    # 绘制主框底部图例
    cursor = x
    for label, fill, stroke in items:
        canvas.circle(cursor, y, 9, fill=fill, stroke=stroke, stroke_width=1.4)
        canvas.text(label, cursor + 18, y + 1, **text_style(15, "500", COLORS["muted"], "start"))
        cursor += 172


def draw_left_panel(canvas: CanvasProtocol) -> None:
    # 绘制左侧 S. cerevisiae essential gene 网络
    x, y, w, h = 90, 190, 510, 560
    draw_panel(canvas, x, y, w, h, "S. cerevisiae essential genes", COLORS["left_header"])
    local_nodes = [
        {"x": 102, "y": 116, "r": 20, "fill": COLORS["essential_blue"], "stroke": "#FFFFFF", "label": "E"},
        {"x": 202, "y": 92, "r": 16, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 308, "y": 122, "r": 19, "fill": COLORS["essential_blue"], "stroke": "#FFFFFF", "label": "E"},
        {"x": 385, "y": 208, "r": 15, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 255, "y": 240, "r": 22, "fill": COLORS["essential_blue"], "stroke": "#FFFFFF", "label": "E"},
        {"x": 126, "y": 270, "r": 15, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 196, "y": 365, "r": 18, "fill": COLORS["essential_blue"], "stroke": "#FFFFFF", "label": "E"},
        {"x": 342, "y": 352, "r": 16, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 420, "y": 324, "r": 18, "fill": COLORS["essential_blue"], "stroke": "#FFFFFF", "label": "E"},
    ]
    edges = [(0, 1), (0, 4), (1, 2), (2, 3), (2, 4), (3, 8), (4, 5), (4, 6), (4, 7), (6, 7), (7, 8)]
    draw_network(canvas, 118, 250, local_nodes, edges)
    draw_legend(canvas, x + 105, y + h - 56, [("Essential", COLORS["essential_blue"], "#FFFFFF"), ("Non-essential", COLORS["nonessential"], COLORS["nonessential_stroke"])])


def draw_gap(canvas: CanvasProtocol) -> None:
    # 绘制中间 transfer / annotation gap 说明区
    x, y, w = 650, 242, 300
    canvas.text("Transfer / annotation gap", x + w / 2, y - 34, **text_style(21, "700", COLORS["ink"], "middle"))
    canvas.line(605, 470, 650, 470, stroke=COLORS["arrow"], stroke_width=2.2, stroke_dasharray="8 8", marker_end="url(#arrowHead)")
    canvas.line(950, 470, 995, 470, stroke=COLORS["arrow"], stroke_width=2.2, stroke_dasharray="8 8", marker_end="url(#arrowHead)")
    for label, box_y in [("uncertain labels", 282), ("incomplete orthology transfer", 407), ("sparse pathogen-specific annotations", 532)]:
        canvas.rect(x, box_y, w, 78, rx=18, ry=18, fill=COLORS["gap_box"], stroke=COLORS["gap_stroke"], stroke_width=1.4)
        canvas.text(label, x + w / 2, box_y + 40, **text_style(17, "600", COLORS["muted"], "middle"))


def draw_right_panel(canvas: CanvasProtocol) -> None:
    # 绘制右侧 Fusarium genome-wide prediction 网络
    x, y, w, h = 1000, 190, 510, 560
    draw_panel(canvas, x, y, w, h, "Fusarium graminearum genome-wide prediction", COLORS["right_header"])
    local_nodes = [
        {"x": 80, "y": 134, "r": 14, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 168, "y": 98, "r": 18, "fill": COLORS["pred_red"], "stroke": "#FFFFFF", "label": "P"},
        {"x": 277, "y": 125, "r": 14, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 390, "y": 110, "r": 20, "fill": COLORS["pred_red"], "stroke": "#FFFFFF", "label": "P"},
        {"x": 115, "y": 260, "r": 16, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 238, "y": 246, "r": 22, "fill": COLORS["pred_red"], "stroke": "#FFFFFF", "label": "P"},
        {"x": 358, "y": 255, "r": 15, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 182, "y": 365, "r": 15, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
        {"x": 305, "y": 373, "r": 19, "fill": COLORS["pred_red"], "stroke": "#FFFFFF", "label": "P"},
        {"x": 420, "y": 350, "r": 14, "fill": COLORS["nonessential"], "stroke": COLORS["nonessential_stroke"]},
    ]
    edges = [(0, 1), (1, 2), (1, 5), (2, 3), (2, 5), (3, 6), (4, 5), (5, 6), (5, 7), (5, 8), (6, 9), (8, 9)]
    draw_network(canvas, 1040, 250, local_nodes, edges)
    draw_legend(canvas, x + 98, y + h - 56, [("Predicted essential", COLORS["pred_red"], "#FFFFFF"), ("Predicted non-essential", COLORS["nonessential"], COLORS["nonessential_stroke"])])


def write_readme() -> None:
    # 写出 Figure 1A 复现说明
    README_OUT.write_text(
        """# Figure 1A Challenge Diagram

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
""",
        encoding="utf-8",
    )


def export_with_rsvg() -> None:
    # 使用 rsvg-convert 从原生 SVG 导出 PDF 和 PNG
    converter = shutil.which("rsvg-convert")
    if converter is None:
        raise RuntimeError("rsvg-convert was not found; cannot export PDF/PNG from SVG.")
    subprocess.run([converter, "-f", "pdf", "-o", str(PDF_OUT), str(SVG_OUT)], check=True)
    subprocess.run([converter, "-f", "png", "-w", "3200", "-h", "1800", "-o", str(PNG_OUT), str(SVG_OUT)], check=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    canvas = make_canvas()
    add_defs(canvas)
    canvas.rect(0, 0, WIDTH, HEIGHT, fill="#FFFFFF")
    draw_title(canvas)
    draw_left_panel(canvas)
    draw_gap(canvas)
    draw_right_panel(canvas)
    canvas.write(SVG_OUT)
    export_with_rsvg()
    write_readme()
    print(f"SVG saved to {SVG_OUT}")
    print(f"PDF saved to {PDF_OUT}")
    print(f"PNG saved to {PNG_OUT}")
    print("Done.")


if __name__ == "__main__":
    main()
