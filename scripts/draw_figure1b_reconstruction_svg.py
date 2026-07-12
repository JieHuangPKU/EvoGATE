#!/usr/bin/env python3
"""
脚本名称: draw_figure1b_reconstruction_svg.py
日期: 2026-05-06
作者: OpenAI/Codex + Jie Huang workflow support
功能描述: 直接生成论文风格 Figure 1B reconstruction 概念示意图原生 SVG，并导出 PDF 和 PNG。
"""

from __future__ import annotations

import html
import shutil
import subprocess
from pathlib import Path
from typing import Protocol


REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "results" / "Figure1B"
SVG_OUT = OUT_DIR / "Figure1B_reconstruction_diagram.svg"
PDF_OUT = OUT_DIR / "Figure1B_reconstruction_diagram.pdf"
PNG_OUT = OUT_DIR / "Figure1B_reconstruction_diagram.png"
README_OUT = OUT_DIR / "README.md"

WIDTH = 1600
HEIGHT = 720

# 统一颜色配置，修改配色优先改这里
COLORS = {
    "ink": "#232323",
    "muted": "#5C6670",
    "faint": "#7A858F",
    "frame": "#BFC9D2",
    "soft_shadow": "#D8DEE7",
    "teal_fill": "#E9F5F4",
    "teal_stroke": "#78AFAE",
    "teal_deep_fill": "#D7ECE8",
    "teal_deep_stroke": "#4F8F8A",
    "final_fill": "#CBE5DF",
    "final_stroke": "#3F7F79",
    "phi_fill": "#FFF2E7",
    "phi_stroke": "#D8A46E",
    "neutral_fill": "#F7F9FA",
    "arrow": "#8F9BA6",
}


class CanvasProtocol(Protocol):
    def add_def(self, text: str) -> None: ...
    def rect(self, x: float, y: float, w: float, h: float, **attrs: object) -> None: ...
    def circle(self, cx: float, cy: float, r: float, **attrs: object) -> None: ...
    def line(self, x1: float, y1: float, x2: float, y2: float, **attrs: object) -> None: ...
    def text(self, text: str, x: float, y: float, **attrs: object) -> None: ...
    def write(self, path: Path) -> None: ...


class SvgCanvas:
    """Small SVG writer used when svgwrite is unavailable."""

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
        attrs = {"x": x, "y": y, "width": w, "height": h, **attrs}
        self.add(f"<rect {format_attrs(attrs)} />")

    def circle(self, cx: float, cy: float, r: float, **attrs: object) -> None:
        attrs = {"cx": cx, "cy": cy, "r": r, **attrs}
        self.add(f"<circle {format_attrs(attrs)} />")

    def line(self, x1: float, y1: float, x2: float, y2: float, **attrs: object) -> None:
        attrs = {"x1": x1, "y1": y1, "x2": x2, "y2": y2, **attrs}
        self.add(f"<line {format_attrs(attrs)} />")

    def text(self, text: str, x: float, y: float, **attrs: object) -> None:
        attrs = {"x": x, "y": y, **attrs}
        self.add(f"<text {format_attrs(attrs)}>{html.escape(text)}</text>")

    def group_start(self, **attrs: object) -> None:
        self.add(f"<g {format_attrs(attrs)}>")

    def group_end(self) -> None:
        self.add("</g>")

    def write(self, path: Path) -> None:
        defs = f"<defs>{''.join(self.defs)}</defs>" if self.defs else ""
        svg = (
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}" role="img" aria-label="Figure 1B reconstruction diagram">\n'
            f"{defs}\n"
            + "\n".join(self.elements)
            + "\n</svg>\n"
        )
        path.write_text(svg, encoding="utf-8")


class SvgwriteCanvas:
    """svgwrite-backed SVG writer that keeps text and shapes editable."""

    def __init__(self, width: int, height: int) -> None:
        import svgwrite

        self.svgwrite = svgwrite
        self.dwg = svgwrite.Drawing(
            filename=str(SVG_OUT),
            size=(width, height),
            viewBox=f"0 0 {width} {height}",
            profile="full",
            debug=False,
        )
        self.dwg.attribs["role"] = "img"
        self.dwg.attribs["aria-label"] = "Figure 1B reconstruction diagram"
        self.arrow_head_iri = "url(#arrowHead)"

    def add_def(self, text: str) -> None:
        # 使用 svgwrite 原生 marker；阴影使用独立浅色偏移图形，不依赖 filter
        marker = self.dwg.marker(
            id="arrowHead",
            insert=(10, 4.5),
            size=(12, 9),
            orient="auto",
            markerUnits="strokeWidth",
        )
        marker.add(self.dwg.path(d="M0,0 L11,4.5 L0,9 z", fill="#B8C0CA"))
        self.dwg.defs.add(marker)

    def rect(self, x: float, y: float, w: float, h: float, **attrs: object) -> None:
        attrs = normalize_svgwrite_attrs(attrs)
        self.dwg.add(self.dwg.rect(insert=(x, y), size=(w, h), **attrs))

    def circle(self, cx: float, cy: float, r: float, **attrs: object) -> None:
        attrs = normalize_svgwrite_attrs(attrs)
        self.dwg.add(self.dwg.circle(center=(cx, cy), r=r, **attrs))

    def line(self, x1: float, y1: float, x2: float, y2: float, **attrs: object) -> None:
        attrs = normalize_svgwrite_attrs(attrs)
        self.dwg.add(self.dwg.line(start=(x1, y1), end=(x2, y2), **attrs))

    def text(self, text: str, x: float, y: float, **attrs: object) -> None:
        attrs = normalize_svgwrite_attrs(attrs)
        self.dwg.add(self.dwg.text(text, insert=(x, y), **attrs))

    def write(self, path: Path) -> None:
        self.dwg.filename = str(path)
        self.dwg.save(pretty=True)


def format_attrs(attrs: dict[str, object]) -> str:
    parts = []
    for key, value in attrs.items():
        if value is None:
            continue
        svg_key = key.replace("_", "-")
        parts.append(f'{svg_key}="{html.escape(str(value), quote=True)}"')
    return " ".join(parts)


def normalize_svgwrite_attrs(attrs: dict[str, object]) -> dict[str, object]:
    converted = {}
    for key, value in attrs.items():
        if value is None:
            continue
        converted[key.replace("_", "-")] = value
    return converted


def make_canvas() -> CanvasProtocol:
    # 优先使用 svgwrite；若当前环境没有该包，则使用本脚本的原生 SVG writer
    try:
        import svgwrite  # noqa: F401

        print("[Figure1B] svgwrite is available; using svgwrite backend.")
        return SvgwriteCanvas(WIDTH, HEIGHT)
    except ImportError:
        print("[Figure1B] svgwrite is not installed; using built-in native SVG writer fallback.")
        return SvgCanvas(WIDTH, HEIGHT)


def add_defs(canvas: CanvasProtocol) -> None:
    # 添加箭头和轻微阴影定义
    canvas.add_def(
        """
        <filter id="softShadow" x="-20%" y="-20%" width="140%" height="140%">
          <feDropShadow dx="0" dy="6" stdDeviation="5" flood-color="#CAD1DB" flood-opacity="0.42"/>
        </filter>
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
    # 绘制左上角 panel label 和左对齐主文 panel 标题
    canvas.text("(B)", 72, 64, **text_style(24, "700", COLORS["ink"], "start"))
    canvas.text(
        "Reconstructing essential-gene labels in a non-model fungal pathogen",
        132,
        64,
        **text_style(22, "600", COLORS["ink"], "start"),
    )
    canvas.text(
        "Cross-species transfer and PHI-base evidence converge into a curated Fusarium positive set",
        132,
        96,
        **text_style(13, "400", COLORS["faint"], "start"),
    )


def draw_node_box(
    canvas: CanvasProtocol,
    x: int,
    y: int,
    w: int,
    h: int,
    title: str,
    subtitle: str,
    fill: str,
    stroke: str,
    weight: str = "500",
) -> None:
    # 绘制统一圆角节点，保持主文 figure 的克制视觉层级
    radius = 12
    canvas.rect(x, y, w, h, rx=radius, ry=radius, fill=fill, stroke=stroke, stroke_width=1.4)
    title_lines = title.split("\n")
    if len(title_lines) == 1:
        title_y_values = [y + 31]
    else:
        title_y_values = [y + 25, y + 43]
    for line, line_y in zip(title_lines, title_y_values):
        canvas.text(line, x + w / 2, line_y, **text_style(15, weight, COLORS["ink"], "middle"))
    if subtitle:
        subtitle_y = y + h - 24
        canvas.text(subtitle, x + w / 2, subtitle_y, **text_style(11, "400", COLORS["faint"], "middle"))


def draw_arrow(canvas: CanvasProtocol, x1: int, y1: int, x2: int, y2: int, dashed: bool = False) -> None:
    # 绘制统一线宽和箭头头部的连接箭头
    attrs: dict[str, object] = {
        "stroke": COLORS["arrow"],
        "stroke_width": 1.5,
        "stroke_linecap": "round",
        "marker_end": "url(#arrowHead)",
    }
    if dashed:
        attrs["stroke_dasharray"] = "5 6"
    canvas.line(x1, y1, x2, y2, **attrs)


def draw_stage_label(canvas: CanvasProtocol, text: str, x: int, y: int) -> None:
    # 绘制轻量级阶段标签，帮助读者区分主链和汇合点
    canvas.text(text, x, y, **text_style(11, "600", COLORS["faint"], "middle"))


def draw_main_chain(canvas: CanvasProtocol) -> None:
    # 绘制 source evidence -> transfer/mapping -> aggregation -> final reconstructed set 主链
    y = 345
    h = 88
    gap = 40
    w = 176
    nodes = [
        (90, y, w, h, "Yeast essential\nsupport", "S. cerevisiae / S. pombe"),
        (90 + (w + gap), y, w, h, "Orthogroup-based\ntransfer", "shared evolutionary context"),
        (90 + 2 * (w + gap), y, w, h, "Canonical Fusarium\ngene mapping", "resolved PH-1 gene IDs"),
        (90 + 3 * (w + gap), y, w, h, "High-confidence\ntransfer positives", "retained transfer evidence"),
    ]
    aggregation = (990, y - 8, 218, 104, "Evidence\naggregation", "transfer + pathogen evidence")
    final = (1276, y - 18, 240, 124, "Reconstructed essential\ngene set", "final curated positives")

    for idx, (x, box_y, box_w, box_h, title, subtitle) in enumerate(nodes):
        draw_node_box(
            canvas,
            x,
            box_y,
            box_w,
            box_h,
            title.replace("\\n", "\n"),
            subtitle,
            COLORS["teal_fill"],
            COLORS["teal_stroke"],
            "500",
        )
        if idx < len(nodes) - 1:
            x2 = nodes[idx + 1][0]
            draw_arrow(canvas, x + box_w + 8, box_y + box_h // 2, x2 - 8, box_y + box_h // 2)

    last_x, last_y, last_w, last_h = nodes[-1][0], nodes[-1][1], nodes[-1][2], nodes[-1][3]
    draw_arrow(canvas, last_x + last_w + 8, last_y + last_h // 2, aggregation[0] - 8, y + h // 2)
    draw_node_box(canvas, *aggregation, fill=COLORS["teal_deep_fill"], stroke=COLORS["teal_deep_stroke"], weight="600")
    draw_arrow(canvas, aggregation[0] + aggregation[2] + 8, y + h // 2, final[0] - 8, y + h // 2)
    draw_node_box(canvas, *final, fill=COLORS["final_fill"], stroke=COLORS["final_stroke"], weight="650")

    # 绘制轻量级阶段标签，形成可读但不喧宾夺主的结构提示
    draw_stage_label(canvas, "source evidence", 184, 296)
    draw_stage_label(canvas, "transfer / mapping", 502, 296)
    draw_stage_label(canvas, "aggregation", 1099, 296)
    draw_stage_label(canvas, "final output", 1396, 296)


def draw_phi_branch(canvas: CanvasProtocol) -> None:
    # 绘制上方并行 PHI-base 分支，汇入 Evidence aggregation
    x, y, w, h = 990, 180, 218, 78
    draw_node_box(
        canvas,
        x,
        y,
        w,
        h,
        "PHI-base lethal-\nsupported genes",
        "pathogen-specific evidence",
        COLORS["phi_fill"],
        COLORS["phi_stroke"],
        "500",
    )
    draw_arrow(canvas, x + w // 2, y + h + 8, x + w // 2, 329, dashed=True)


def draw_context_notes(canvas: CanvasProtocol) -> None:
    # 绘制底部简短说明，避免堆信息但保留科学语义
    canvas.text(
        "Cross-species essentiality evidence is filtered through orthogroups and canonical gene mapping, then merged with pathogen-specific lethal evidence.",
        90,
        560,
        **text_style(12, "400", COLORS["muted"], "start"),
    )
    canvas.text(
        "Goal: a compact reconstructed positive set for genome-wide essential-gene prediction in Fusarium graminearum.",
        90,
        587,
        **text_style(12, "400", COLORS["muted"], "start"),
    )


def write_readme() -> None:
    # 写出 Figure 1B 复现说明
    README_OUT.write_text(
        f"""# Figure 1B Reconstruction Diagram

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
""",
        encoding="utf-8",
    )


def export_with_rsvg() -> None:
    # 使用 rsvg-convert 从原生 SVG 导出 PDF 和 PNG
    converter = shutil.which("rsvg-convert")
    if converter is None:
        raise RuntimeError("rsvg-convert was not found; cannot export PDF/PNG from SVG.")

    subprocess.run([converter, "-f", "pdf", "-o", str(PDF_OUT), str(SVG_OUT)], check=True)
    subprocess.run([converter, "-f", "png", "-w", "3200", "-h", "1440", "-o", str(PNG_OUT), str(SVG_OUT)], check=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # 生成原生 SVG 画布
    canvas = make_canvas()
    add_defs(canvas)
    canvas.rect(0, 0, WIDTH, HEIGHT, fill="#FFFFFF")

    # 绘制 Figure 1B 各组成部分
    draw_title(canvas)
    draw_phi_branch(canvas)
    draw_main_chain(canvas)
    draw_context_notes(canvas)

    # 输出 SVG、PDF、PNG 和 README
    canvas.write(SVG_OUT)
    export_with_rsvg()
    write_readme()

    print(f"SVG saved to {SVG_OUT}")
    print(f"PDF saved to {PDF_OUT}")
    print(f"PNG saved to {PNG_OUT}")
    print("Done.")


if __name__ == "__main__":
    main()
