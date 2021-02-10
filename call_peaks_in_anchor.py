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
import dask.array as da
import tifffile as tf
import trackpy as tp
import pysnooper


@pysnooper.snoop()
def main(args):
    # Call peaks
    processed_anchor = da.from_zarr(args.anchor_zarr, "processed")
    # processed_anchor = np.transpose(processed_anchor, (2, 0, 1))

    if processed_anchor.shape[0] > 1:
        tracks = tp.link_df(
                tp.batch(np.array(processed_anchor),
                    diameter=args.spot_diameter,
                    percentile=args.trackpy_percentile,
                    separation=args.peak_separation
                    ),
                search_range=5
            )
        # peaks = tracks['particle'][tracks['frame'] == 0].unique()
    else:
        peaks = tp.locate(
            np.array(processed_anchor),
            diameter=args.spot_diameter,
            percentile=args.trackpy_percentile,
            separation=args.peak_separation
        )

    peaks.to_csv("%s_peaks.tsv" %args.stem, sep="\t", index=False)


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
