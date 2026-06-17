#!/usr/bin/env python3
"""Crop an 8-bit non-interlaced PNG to a vertical row range.

This tiny utility is included because the slide 13b pedigree source has no
local regeneration recipe and the environment does not provide ImageMagick.
It intentionally implements only the PNG subset used by the deck assets:
8-bit RGB/RGBA, non-interlaced PNGs.
"""

from __future__ import annotations

import argparse
import binascii
import struct
import zlib
from pathlib import Path


PNG_SIG = b"\x89PNG\r\n\x1a\n"


def paeth(a: int, b: int, c: int) -> int:
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def read_chunks(blob: bytes):
    if not blob.startswith(PNG_SIG):
        raise ValueError("not a PNG file")
    offset = len(PNG_SIG)
    while offset < len(blob):
        if offset + 8 > len(blob):
            raise ValueError("truncated PNG chunk")
        length = struct.unpack(">I", blob[offset:offset + 4])[0]
        ctype = blob[offset + 4:offset + 8]
        start = offset + 8
        end = start + length
        data = blob[start:end]
        yield ctype, data
        offset = end + 4
        if ctype == b"IEND":
            break


def unfilter_scanlines(raw: bytes, width: int, height: int, channels: int) -> list[bytes]:
    bpp = channels
    stride = width * channels
    rows: list[bytes] = []
    offset = 0
    prev = bytearray(stride)
    for _ in range(height):
        filter_type = raw[offset]
        offset += 1
        scan = bytearray(raw[offset:offset + stride])
        offset += stride
        out = bytearray(stride)
        for i, value in enumerate(scan):
            left = out[i - bpp] if i >= bpp else 0
            up = prev[i]
            upper_left = prev[i - bpp] if i >= bpp else 0
            if filter_type == 0:
                recon = value
            elif filter_type == 1:
                recon = value + left
            elif filter_type == 2:
                recon = value + up
            elif filter_type == 3:
                recon = value + ((left + up) // 2)
            elif filter_type == 4:
                recon = value + paeth(left, up, upper_left)
            else:
                raise ValueError(f"unsupported PNG filter type {filter_type}")
            out[i] = recon & 0xFF
        rows.append(bytes(out))
        prev = out
    return rows


def chunk(ctype: bytes, data: bytes) -> bytes:
    crc = binascii.crc32(ctype)
    crc = binascii.crc32(data, crc) & 0xFFFFFFFF
    return struct.pack(">I", len(data)) + ctype + data + struct.pack(">I", crc)


def crop_region(input_path: Path, output_path: Path, y_offset: int, height: int) -> None:
    blob = input_path.read_bytes()
    ihdr = None
    idat_parts = []
    for ctype, data in read_chunks(blob):
        if ctype == b"IHDR":
            ihdr = data
        elif ctype == b"IDAT":
            idat_parts.append(data)
    if ihdr is None:
        raise ValueError("PNG is missing IHDR")

    width, old_height, bit_depth, color_type, comp, filt, interlace = struct.unpack(
        ">IIBBBBB", ihdr
    )
    if bit_depth != 8 or color_type not in (2, 6) or comp != 0 or filt != 0 or interlace != 0:
        raise ValueError(
            "unsupported PNG format; expected 8-bit non-interlaced RGB/RGBA"
        )
    if y_offset < 0:
        raise ValueError("y-offset must be non-negative")
    if not 1 <= height <= old_height - y_offset:
        raise ValueError(f"height must be between 1 and {old_height - y_offset}")
    channels = 3 if color_type == 2 else 4
    rows = unfilter_scanlines(zlib.decompress(b"".join(idat_parts)),
                              width, old_height, channels)
    cropped_rows = rows[y_offset:(y_offset + height)]
    cropped = b"".join(b"\x00" + row for row in cropped_rows)
    new_ihdr = struct.pack(
        ">IIBBBBB", width, height, bit_depth, color_type, comp, filt, interlace
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(
        PNG_SIG +
        chunk(b"IHDR", new_ihdr) +
        chunk(b"IDAT", zlib.compress(cropped, level=9)) +
        chunk(b"IEND", b"")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--height", type=int, required=True)
    parser.add_argument("--y-offset", type=int, default=0)
    args = parser.parse_args()
    crop_region(args.input, args.output, args.y_offset, args.height)


if __name__ == "__main__":
    main()
