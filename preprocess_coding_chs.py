#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Preprocess coding channel images
"""
import argparse
import dask.array as da
import numpy as np
from skimage.restoration import denoise_wavelet
from skimage.morphology import white_tophat, disk
import zarr


def main(args):
    imgs = da.from_zarr(args.zarr, "coding_ch_images")

    # white tophat enhancement
    hat_enhenced = imgs.map_overlap(
        white_tophat,
        dtype=np.uint16,
        depth=(args.overlaps, args.overlaps, 0, 0),
        selem=np.expand_dims(disk(args.spot_diameter), (-2, -1)),
    )

    # Wavelet denoise
    denoised = (
        hat_enhenced.map_overlap(
            denoise_wavelet,
            depth=(args.overlaps, args.overlaps, 0, 0),
            method="BayesShrink",
            mode="soft",
            sigma=args.spot_diameter * 1.5,
            rescale_sigma=True,
            multichannel=False,
        )
        * 10 ** 4
    ).astype(np.uint16)

    rechunked = denoised.rechunk({0: -1, 1: -1, 2: 1, 3: 1})
    rechunked = rechunked.rechunk({0: "auto", 1: "auto", 2: 1, 3: 1})
    rechunked.to_zarr(args.out, "normalized_coding_chs")

    source = zarr.open(args.zarr, mode="r")
    store = zarr.DirectoryStore(args.out)
    g = zarr.group(store=store)
    zarr.copy(source["cycle"], g)
    zarr.copy(source["base"], g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-zarr", type=str, required=True)
    parser.add_argument("-out", type=str, required=True)
    parser.add_argument("-spot_diameter", type=int, default=5)
    # parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-quantile_for_norm", type=float, default=99.99)

    args = parser.parse_args()

    main(args)
