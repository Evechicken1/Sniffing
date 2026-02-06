# -*- coding: utf-8 -*-
"""
Script to merge TDMS files where the time is reset between files

@author: Xander
"""

import os
import glob
import numpy as np
import pandas as pd
from nptdms import TdmsFile, TdmsWriter, ChannelObject

# --- CONFIG ---
INPUT_DIR = r"Z:\Users\Xander\TDMS_test\25_08_18_XT002"   # folder containing the .tdms files to append
OUTPUT_TDMS = os.path.join(INPUT_DIR, "merged.tdms")
WRITE_CSV = True                               # set True to also export CSVs with explicit time columns
CSV_DIR = os.path.join(INPUT_DIR, "merged_csv") # where CSVs go if WRITE_CSV=True
# --------------

def natural_key(s):
    # Natural sort helper: "file2" < "file10"
    import re
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]

def list_tdms_files(input_dir):
    files = glob.glob(os.path.join(input_dir, "*.tdms"))
    if not files:
        raise FileNotFoundError(f"No .tdms files found in: {input_dir}")
    files.sort(key=natural_key)
    return files

def get_structure_and_props(tdms_file):
    """Return {group: {channel: {'wf_increment':..., 'wf_start_time':..., 'unit_string':...}}}"""
    info = {}
    tf = TdmsFile.read(tdms_file)
    for g in tf.groups():
        gname = g.name
        info[gname] = {}
        for ch in g.channels():
            cname = ch.name
            p = ch.properties
            info[gname][cname] = {
                "wf_increment": p.get("wf_increment", None),
                "wf_start_time": p.get("wf_start_time", None),
                "unit_string": p.get("unit_string", None),
            }
    return info

def check_same_structure(struct_list):
    """Ensure all files have the same groups/channels; return canonical (from first)."""
    ref = struct_list[0]
    for i, s in enumerate(struct_list[1:], start=2):
        if set(s.keys()) != set(ref.keys()):
            raise ValueError(f"Group mismatch between file1 and file{i}.")
        for g in ref:
            if set(s[g].keys()) != set(ref[g].keys()):
                raise ValueError(f"Channel mismatch in group '{g}' between file1 and file{i}.")
    return ref

def check_sampling(struct_list):
    """Verify wf_increment is consistent for each channel across files."""
    ref = struct_list[0]
    for i, s in enumerate(struct_list[1:], start=2):
        for g in ref:
            for c in ref[g]:
                inc0 = ref[g][c]["wf_increment"]
                incN = s[g][c]["wf_increment"]
                if (inc0 is not None) and (incN is not None) and (inc0 != incN):
                    raise ValueError(
                        f"Sampling mismatch for {g}/{c}: file1 wf_increment={inc0}, file{i} wf_increment={incN}"
                    )

def merge_tdms_folder(input_dir, output_tdms, write_csv=False, csv_dir=None):
    files = list_tdms_files(input_dir)

    # Inspect structures and properties
    structs = [get_structure_and_props(f) for f in files]
    canon = check_same_structure(structs)
    check_sampling(structs)

    if write_csv and csv_dir:
        os.makedirs(csv_dir, exist_ok=True)

    # Read & concatenate channel data across files
    # We’ll write out in segments per group to the resulting TDMS
    with TdmsWriter(output_tdms) as writer:
        for gname in canon.keys():
            # We’ll build ChannelObjects for all channels in this group, then write in one segment
            group_channel_objs = []

            for cname, p0 in canon[gname].items():
                data_pieces = []
                # collect data from each file
                for path in files:
                    t = TdmsFile.read(path)
                    data_pieces.append(t[gname][cname].data)

                merged = np.concatenate(data_pieces)

                # Reuse reference properties (wf_start_time of first file makes a continuous implicit time axis)
                props = {}
                if p0["wf_start_time"] is not None:
                    props["wf_start_time"] = p0["wf_start_time"]
                if p0["wf_increment"] is not None:
                    props["wf_increment"] = p0["wf_increment"]
                if p0["unit_string"] is not None:
                    props["unit_string"] = p0["unit_string"]

                # Make the channel object
                ch_obj = ChannelObject(gname, cname, merged, properties=props)
                group_channel_objs.append(ch_obj)

                # Optional CSV export with explicit time column
                if write_csv:
                    # If wf_increment known, build continuous time from start at 0 with the increment.
                    # If not, try to derive from the first file's time_track() (if available), else just index.
                    dt = p0["wf_increment"]
                    if dt is not None:
                        # NI stores seconds as a float; we’ll export time in seconds
                        time = np.arange(merged.size, dtype=float) * float(dt)
                    else:
                        # Try to reconstruct from first file’s channel time track
                        try:
                            first_ch = TdmsFile.read(files[0])[gname][cname]
                            tt = getattr(first_ch, "time_track", None)
                            if callable(tt):
                                # Build per-file times, then offset each subsequent file’s time
                                time_parts = []
                                offset = 0.0
                                for path in files:
                                    ch = TdmsFile.read(path)[gname][cname]
                                    ttrack = ch.time_track()  # numpy datetime64
                                    # Convert to seconds starting at 0 within the file
                                    t0 = ttrack[0]
                                    secs = (ttrack - t0).astype("timedelta64[ns]").astype(np.int64) / 1e9
                                    # offset by cumulative end
                                    secs = secs + offset
                                    offset = secs[-1] + (secs[1] - secs[0] if secs.size > 1 else 0.0)
                                    time_parts.append(secs)
                                time = np.concatenate(time_parts)
                            else:
                                time = np.arange(merged.size, dtype=float)
                        except Exception:
                            time = np.arange(merged.size, dtype=float)

                    df = pd.DataFrame({"time_s": time, cname: merged})
                    out_csv = os.path.join(csv_dir, f"{gname}__{cname}.csv")
                    df.to_csv(out_csv, index=False)

            # write the whole group as a segment (efficient)
            writer.write_segment(group_channel_objs)

    return output_tdms

if __name__ == "__main__":
    out = merge_tdms_folder(INPUT_DIR, OUTPUT_TDMS, write_csv=WRITE_CSV, csv_dir=CSV_DIR if WRITE_CSV else None)
    print(f"✅ Merged TDMS written to: {out}")
    if WRITE_CSV:
        print(f"✅ CSVs written to: {CSV_DIR}")

