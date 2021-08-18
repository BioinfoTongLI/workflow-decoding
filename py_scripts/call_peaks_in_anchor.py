#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Preprocess all channels for decoding
"""
import argparse
import numpy as np
import trackpy as tp
import pandas as pd
import xarray as xr
from concurrent.futures import ThreadPoolExecutor


def trackpy_locate(ch_name, img):
    peaks = tp.locate(
        np.array(img),
        diameter=args.spot_diameter,
        percentile=args.trackpy_percentile,
        separation=args.peak_separation,
    )
    peaks["ch_name"] = ch_name
    peaks["frame"] = ch_name[2]
    return peaks


def main(args):
    processed_anchor = xr.open_zarr(args.anchor_zarr)

    # print(processed_anchor.data_vars)
    print(len(processed_anchor.data_vars))

    # print(processed_anchor.to_array().data.rechunk((-1, "auto", "auto")))
    # rechunked = processed_anchor.to_array().data.rechunk((1, -1, -1))
    # print(rechunked.map_blocks(trackpy_locate))
    # print(xr.map_blocks(trackpy_locate, ))
    peak_calling_futures = []
    # Call peaks
    for v in processed_anchor.data_vars:
        print(v, processed_anchor.data_vars[v])
        executor = ThreadPoolExecutor(15)
        peak_calling_futures.append(
            executor.submit(trackpy_locate, v, processed_anchor.data_vars[v])
            # trackpy_locate(v, processed_anchor.data_vars[v])
        )
    # print(peak_calling_futures)
    peaks = [obj.result() for obj in peak_calling_futures]
    concat_peaks = pd.concat(peaks)
    tracked_df = tp.link(concat_peaks, search_range=5)
    print(tracked_df)
    tracked_df.to_csv("%s_peaks.tsv" % args.stem, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-anchor_zarr", type=str, required=True)
    parser.add_argument("-stem", type=str, required=True)
    parser.add_argument("-spot_diameter", type=int, default=5)
    parser.add_argument("-trackpy_percentile", type=int, default=64)
    parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-trackpy_search_range", type=int, default=9)
    parser.add_argument("-peak_separation", type=int, default=5)

    args = parser.parse_args()

    main(args)
