#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Preprocess cdoing channel images
"""
import argparse
import dask.array as da
import pysnooper
import numpy as np
from skimage.restoration import denoise_wavelet
from skimage.morphology import white_tophat, disk
import zarr


def normalize(stack, quantile):
    img_min = np.nanmin(stack, axis=(0, 1), keepdims=True)
    img_max = np.nanpercentile(stack, q=quantile, axis=(0, 1), keepdims=True)
    normalized = (stack - img_min) / (img_max - img_min)
    return normalized


def normalize_gpu(stack, quantile):
    import cupy

    stack = cupy.array(np.array(stack))
    img_min = cupy.percentile(stack, q=quantile, axis=(0, 1), keepdims=True)
    img_max = cupy.percentile(stack, q=1 - quantile, axis=(0, 1), keepdims=True)
    return (stack - img_min) / (img_max - img_min)


@pysnooper.snoop()
def main(args):
    imgs = da.from_zarr(args.zarr, "coding_ch_images")

    # white tophat enhancement
    hat_enhenced = imgs.map_overlap(
        white_tophat,
        dtype=np.uint16,
        depth=(args.overlaps, args.overlaps, 0, 0),
        selem=np.expand_dims(disk(args.whitehat_disk_diam), (-2, -1)),
    )

    # Wavelet denoise
    denoised = hat_enhenced.map_overlap(
            denoise_wavelet,
            depth=(args.overlaps, args.overlaps, 0, 0),
            method="BayesShrink",
            mode="soft",
            sigma=5,
            rescale_sigma=True,
            multichannel=False,
        ) * 10 ** 4
    denoised = denoised.astype(np.uint16)

    # nomralization
    # rechunked = imgs.rechunk(
    rechunked = denoised.rechunk({0: -1, 1: -1, 2: 1, 3: 1})
    normalized = rechunked.map_blocks(normalize, quantile=args.quantile_for_norm)
    # normalized = rechunked.map_blocks(normalize_gpu, args.quantile_for_norm, dtype="float32")
    rechunked_normalized = normalized.rechunk({0: "auto", 1: "auto", 2: 1, 3: 1})
    rechunked_normalized.to_zarr(args.out, "normalized_coding_chs")

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
    parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-quantile_for_norm", type=float, default=99.99)

    args = parser.parse_args()

    main(args)
