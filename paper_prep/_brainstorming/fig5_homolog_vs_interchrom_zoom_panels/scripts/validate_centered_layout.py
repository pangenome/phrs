#!/usr/bin/env python3
"""Validate that the rendered Fig5 zoom panel stack is centered on the page."""

from __future__ import annotations

import struct
import sys
import zlib
from pathlib import Path


PNG_PATH = Path(
    "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/"
    "fig5_homolog_vs_interchrom_zoom_panels.png"
)
MAX_CENTER_OFFSET_PX = 20
WHITE_THRESHOLD = 248


def paeth_predictor(a: int, b: int, c: int) -> int:
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def read_png_rgb(path: Path) -> tuple[int, int, list[bytes]]:
    with path.open("rb") as handle:
        if handle.read(8) != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"{path} is not a PNG")

        width = height = bit_depth = color_type = None
        compressed = bytearray()
        while True:
            length_bytes = handle.read(4)
            if not length_bytes:
                raise ValueError(f"{path} ended before IEND")
            length = struct.unpack(">I", length_bytes)[0]
            chunk_type = handle.read(4)
            payload = handle.read(length)
            handle.read(4)

            if chunk_type == b"IHDR":
                width, height, bit_depth, color_type, _comp, _filter, interlace = struct.unpack(
                    ">IIBBBBB", payload
                )
                if bit_depth != 8 or color_type != 2 or interlace != 0:
                    raise ValueError("validator expects non-interlaced 8-bit RGB PNG output")
            elif chunk_type == b"IDAT":
                compressed.extend(payload)
            elif chunk_type == b"IEND":
                break

        if width is None or height is None or bit_depth is None or color_type is None:
            raise ValueError(f"{path} is missing IHDR")

    raw = zlib.decompress(bytes(compressed))
    bytes_per_pixel = 3
    stride = width * bytes_per_pixel
    rows: list[bytes] = []
    previous = bytearray(stride)
    index = 0

    for _y in range(height):
        filter_type = raw[index]
        index += 1
        scanline = raw[index : index + stride]
        index += stride
        reconstructed = bytearray(stride)

        for x in range(stride):
            left = reconstructed[x - bytes_per_pixel] if x >= bytes_per_pixel else 0
            up = previous[x]
            up_left = previous[x - bytes_per_pixel] if x >= bytes_per_pixel else 0

            if filter_type == 0:
                value = scanline[x]
            elif filter_type == 1:
                value = scanline[x] + left
            elif filter_type == 2:
                value = scanline[x] + up
            elif filter_type == 3:
                value = scanline[x] + ((left + up) // 2)
            elif filter_type == 4:
                value = scanline[x] + paeth_predictor(left, up, up_left)
            else:
                raise ValueError(f"unsupported PNG filter type {filter_type}")

            reconstructed[x] = value & 0xFF

        rows.append(bytes(reconstructed))
        previous = reconstructed

    return width, height, rows


def main() -> int:
    width, _height, rows = read_png_rgb(PNG_PATH)
    nonwhite_x: list[int] = []

    for row in rows:
        for x in range(width):
            r, g, b = row[x * 3 : x * 3 + 3]
            if r <= WHITE_THRESHOLD or g <= WHITE_THRESHOLD or b <= WHITE_THRESHOLD:
                nonwhite_x.append(x)

    if not nonwhite_x:
        print(f"{PNG_PATH}: no visible content found", file=sys.stderr)
        return 1

    left = min(nonwhite_x)
    right = max(nonwhite_x)
    content_center = (left + right) / 2
    canvas_center = (width - 1) / 2
    offset = abs(content_center - canvas_center)

    print(
        f"{PNG_PATH}: width={width} left={left} right={right} "
        f"content_center={content_center:.1f} canvas_center={canvas_center:.1f} "
        f"offset={offset:.1f}px"
    )

    if offset > MAX_CENTER_OFFSET_PX:
        print(
            f"center offset {offset:.1f}px exceeds {MAX_CENTER_OFFSET_PX}px tolerance",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
