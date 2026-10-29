#!/usr/bin/env python3
"""Writes the LocusMasker logo and favicons in this folder from one definition (numpy only; no browser or image library).

    python viewer/icons/make_icons.py

A 3 x 3 grid of loci by sequences, three of them masked, laid on the Lo Shu square (4 9 2 / 3 5 7 / 8 1 6): red (too
divergent) on 1, blue (copied) on 2, yellow (both) on 3, the values the report's side bits give those outcomes.

    logo.svg        the logo at any size: cells 10, gaps 1, no margin, on a 32-unit square
    favicon-32.png  the same at 32 px, every edge on a whole pixel (what a high-density tab bar draws)
    favicon-16.svg  the 16 px drawing: cells 10, 8, 10 with gaps 2, so at 16 px cells of 5, 4, 5 pixels and 1-pixel gaps
    favicon-16.png  that drawing at 16 px
"""
import os, struct, zlib, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GREY, RED, BLUE, YELLOW = "#c4c2bd", "#b5172e", "#2f6bd8", "#e9ca00"
MASKED = {(2, 1): RED, (0, 2): BLUE, (1, 0): YELLOW}        # (row, column) of Lo Shu 1, 2 and 3
RX = 1.0                                                    # corner radius, in units of the 32-unit square
DRAWINGS = {"logo": ([10, 10, 10], 1), "favicon-16": ([10, 8, 10], 2)}   # cell sizes along each axis, gap; no margin

def cells(sizes, gap):
    at = np.concatenate([[0], np.cumsum(np.array(sizes[:-1]) + gap)])
    return [(at[k], at[r], sizes[k], sizes[r], MASKED.get((r, k), GREY)) for r in range(3) for k in range(3)]

def svg(sizes, gap):
    rects = "".join(f'<rect x="{x:g}" y="{y:g}" width="{w:g}" height="{h:g}" rx="{RX:g}" fill="{c}"/>' for x, y, w, h, c in cells(sizes, gap))
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32">{rects}</svg>\n'

def raster(sizes, gap, n, ss=16):
    """RGBA pixels at n x n: each pixel's coverage of each rounded cell from ss x ss samples"""
    u = (np.arange(n * ss) + 0.5) / (n * ss) * 32.0                 # sample positions in units
    X, Y = np.meshgrid(u, u)
    rgba = np.zeros((n, n, 4))
    for x, y, w, h, c in cells(sizes, gap):
        cx = np.clip(X, x + RX, x + w - RX); cy = np.clip(Y, y + RX, y + h - RX)
        inside = (X >= x) & (X <= x + w) & (Y >= y) & (Y <= y + h) & ((X - cx) ** 2 + (Y - cy) ** 2 <= RX * RX)
        cov = inside.reshape(n, ss, n, ss).mean(axis=(1, 3))
        rgb = np.array([int(c[i:i + 2], 16) for i in (1, 3, 5)], float)
        rgba[..., :3] += cov[..., None] * rgb; rgba[..., 3] += cov  # cells never share a pixel: the gaps are at least 1 px
    a = rgba[..., 3]; rgb = np.where(a[..., None] > 0, rgba[..., :3] / np.maximum(a[..., None], 1e-9), 0)
    return np.dstack([rgb, a * 255]).round().clip(0, 255).astype(np.uint8)

def png(px):
    h, w, _ = px.shape
    raw = b"".join(b"\x00" + px[r].tobytes() for r in range(h))
    chunk = lambda t, d: struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b"")

if __name__ == "__main__":
    for name, (sizes, gap) in DRAWINGS.items():
        open(os.path.join(HERE, name + ".svg"), "w").write(svg(sizes, gap))
    open(os.path.join(HERE, "favicon-32.png"), "wb").write(png(raster(*DRAWINGS["logo"], 32)))
    open(os.path.join(HERE, "favicon-16.png"), "wb").write(png(raster(*DRAWINGS["favicon-16"], 16)))
    print("wrote", ", ".join(sorted(f for f in os.listdir(HERE) if f != os.path.basename(__file__))))
