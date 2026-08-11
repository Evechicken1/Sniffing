"""Convert a FLIR ResearchIR .seq file to a grayscale mp4.

FLIR .seq files recorded through ResearchIR are a concatenation of FFF frame
records. Each frame carries an index table; the raw 16-bit thermal image lives
in the record whose mainType == 1. This script parses those records directly
(no FLIR SDK / exiftool needed), normalises the 16-bit counts to 8-bit
grayscale, and writes an mp4.

Usage:
    python seq_to_mp4.py INPUT.seq [-o OUTPUT.mp4] [--fps 60]
                         [--per-frame] [--pmin 0.5] [--pmax 99.5]

Notes:
- Default normalisation is *global* (one contrast window across the whole clip,
  taken from the --pmin/--pmax percentiles of a frame sample). This preserves
  relative temperature changes and avoids per-frame brightness flicker. Use
  --per-frame to maximise contrast on each frame independently instead.
- The .seq has no readable per-frame timestamp, so --fps only sets playback
  speed; set it to the camera's acquisition rate for real-time playback.
"""

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import imageio.v2 as imageio

FFF_MAGIC = b"FFF\x00"
RAW_IMAGE_MAINTYPE = 0x0001  # index record holding the raw thermal image


def iter_frames(data):
    """Yield (width, height, uint16_image) for each FFF frame in the buffer."""
    base = 0
    n = len(data)
    while base < n and data[base : base + 4] == FFF_MAGIC:
        index_off = struct.unpack_from("<I", data, base + 0x18)[0]
        n_records = struct.unpack_from("<I", data, base + 0x1C)[0]

        frame_end = 0
        image_rec = None
        for k in range(n_records):
            o = base + index_off + k * 32
            main_type = struct.unpack_from("<H", data, o)[0]
            data_off, data_size = struct.unpack_from("<II", data, o + 12)
            frame_end = max(frame_end, data_off + data_size)
            if main_type == RAW_IMAGE_MAINTYPE:
                image_rec = (base + data_off, data_size)

        if image_rec is None:
            raise ValueError(f"No raw-image record in FFF frame at offset {base}")

        rec_off, _ = image_rec
        # 32-byte image sub-header: width @ +2, height @ +4 (uint16 LE), 16bpp.
        width = struct.unpack_from("<H", data, rec_off + 2)[0]
        height = struct.unpack_from("<H", data, rec_off + 4)[0]
        pixels_off = rec_off + 32
        count = width * height
        img = np.frombuffer(data, dtype="<u2", count=count, offset=pixels_off)
        yield width, height, img.reshape(height, width)

        base += frame_end


def load_frames(path):
    """Read all frames from a .seq file into a (n, h, w) uint16 array."""
    data = Path(path).read_bytes()
    frames = []
    dims = None
    for w, h, img in iter_frames(data):
        if dims is None:
            dims = (w, h)
        elif (w, h) != dims:
            raise ValueError(f"Frame size changed within file: {dims} -> {(w, h)}")
        frames.append(img)
    if not frames:
        raise ValueError(f"No FFF frames found in {path} (not a ResearchIR .seq?)")
    return np.stack(frames)


def to_uint8(stack, per_frame, pmin, pmax):
    """Normalise a uint16 stack to uint8 grayscale."""
    stack = stack.astype(np.float32)
    if per_frame:
        out = np.empty(stack.shape, dtype=np.uint8)
        for i, frame in enumerate(stack):
            lo, hi = np.percentile(frame, [pmin, pmax])
            out[i] = _scale(frame, lo, hi)
        return out
    # Global window from a sample of frames (fast, robust to outliers).
    sample = stack[:: max(1, len(stack) // 50)]
    lo, hi = np.percentile(sample, [pmin, pmax])
    return _scale(stack, lo, hi)


def _scale(arr, lo, hi):
    if hi <= lo:
        hi = lo + 1.0
    scaled = (arr - lo) / (hi - lo)
    return np.clip(scaled * 255.0, 0, 255).astype(np.uint8)


def convert(input_path, output_path=None, fps=60, per_frame=False,
            pmin=0.5, pmax=99.5):
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_suffix(".mp4")
    else:
        output_path = Path(output_path)
        # Treat an extensionless / non-video path as a filename and add .mp4.
        if output_path.suffix.lower() not in {".mp4", ".mov", ".avi", ".mkv"}:
            output_path = output_path.with_name(output_path.name + ".mp4")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading {input_path.name} ...", flush=True)
    stack = load_frames(input_path)
    n, h, w = stack.shape
    print(f"  {n} frames, {w}x{h}, 16-bit")

    print("  normalising to 8-bit grayscale "
          f"({'per-frame' if per_frame else 'global'} {pmin}-{pmax} pct) ...",
          flush=True)
    gray = to_uint8(stack, per_frame, pmin, pmax)

    print(f"  encoding -> {output_path} @ {fps} fps ...", flush=True)
    # macro_block_size=8 lets odd 320x240-style sizes through without padding.
    with imageio.get_writer(output_path, fps=fps, codec="libx264",
                            quality=8, macro_block_size=8) as writer:
        for frame in gray:
            writer.append_data(frame)
    print(f"Done: {output_path}")
    return output_path


def main(argv=None):
    p = argparse.ArgumentParser(description="Convert a FLIR ResearchIR .seq to mp4.")
    p.add_argument("input", help="Path to the .seq file")
    p.add_argument("-o", "--output", help="Output .mp4 path (default: alongside input)")
    p.add_argument("--fps", type=float, default=60, help="Playback frame rate (default 60)")
    p.add_argument("--per-frame", action="store_true",
                   help="Normalise each frame independently (max contrast, may flicker)")
    p.add_argument("--pmin", type=float, default=0.5, help="Lower percentile for contrast")
    p.add_argument("--pmax", type=float, default=99.5, help="Upper percentile for contrast")
    args = p.parse_args(argv)
    convert(args.input, args.output, args.fps, args.per_frame, args.pmin, args.pmax)


if __name__ == "__main__":
    sys.exit(main())
