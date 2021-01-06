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
from skimage.restoration import denoise_wavelet
from skimage.morphology import white_tophat, disk
import dask.array as da
from dask_image import imread
import tifffile as tf
from apeer_ometiff_library import omexmlClass
import trackpy as tp


def get_processed_anchor_ch(anchor_img, chunk_size):
    anchor_img = da.transpose(anchor_img, (1, 2, 0))
    anchor_img = anchor_img.rechunk((chunk_size, chunk_size, 1))

    # Tophat enhancement
    hat_enhenced = anchor_img.map_blocks(
        white_tophat,
        selem=np.expand_dims(disk(args.whitehat_disk_diam), -1),
        dtype=np.uint16,
    )

    # Wavelet denoise
    denoised = hat_enhenced.map_overlap(
        denoise_wavelet,
        depth=(args.overlaps, args.overlaps, 0),
        trim=True,
        boundary="nearest",
        method="BayesShrink",
        mode="soft",
        sigma=5,
        rescale_sigma=True,
        multichannel=False,
    ) * 10 ** 4

    return denoised.astype(np.uint16).compute()


def main(args):
    # Retrieve anchor channel
    imgs = imread.imread(args.ome_tif)
    with tf.TiffFile(args.ome_tif) as fh:
        metadata = omexmlClass.OMEXML(fh.ome_metadata)
    pixels = metadata.image(0).Pixels
    ch_names = np.array(pixels.get_channel_names())

    processed_anchor = get_processed_anchor_ch(imgs[ch_names == args.known_anchor],
            args.chunk_size)

    # Save
    da.array(processed_anchor).to_zarr("anchor_image.zarr", "processed_anchor")
    # Call peaks
    peaks = tp.locate(
        processed_anchor,
        diameter=args.spot_diameter,
        percentile=args.trackpy_percentile,
    )

    peaks.to_csv("anchor_peaks.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-ome_tif", type=str, required=True)
    parser.add_argument("-spot_diameter", type=int, default=5)
    parser.add_argument("-trackpy_percentile", type=int, default=64)
    parser.add_argument("-known_anchor", type=str, default="c01 anchor")
    parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-chunk_size", type=int, default=4000)

    args = parser.parse_args()

    main(args)
