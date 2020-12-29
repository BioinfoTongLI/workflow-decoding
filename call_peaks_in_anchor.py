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


def get_processed_anchor_ch(anchor_img):
    chunk_size = args.chunk_size
    anchor_img = da.transpose(anchor_img, (1, 2, 0))
    anchor_img = anchor_img.rechunk((chunk_size, chunk_size, 1))

    # Tophat enhancement
    hat_enhenced = da.map_blocks(
        white_tophat,
        anchor_img,
        selem=np.expand_dims(disk(args.whitehat_disk_diam), -1),
        dtype=np.uint16,
    )

    # Wavelet denoise
    denoised = da.map_overlap(
        denoise_wavelet,
        hat_enhenced,
        depth=(args.overlaps, args.overlaps, 0),
        trim=True,
        boundary="nearest",
        method="BayesShrink",
        mode="soft",
        sigma=151,
        rescale_sigma=True,
        multichannel=False,
    )
    # Save
    # tf.imwrite("processed_anchor.tif", (denoised.compute() * 10 ** 4).astype(np.uint16))
    return (denoised.compute() * 10 ** 4).astype(np.uint16)


def main(args):
    # Retrieve anchor channel
    imgs = imread.imread(args.ome_tif)
    with tf.TiffFile(args.ome_tif) as fh:
        metadata = omexmlClass.OMEXML(fh.ome_metadata)
    pixels = metadata.image(0).Pixels
    ch_names = np.array(pixels.get_channel_names())

    processed_anchor = get_processed_anchor_ch(imgs[ch_names == args.known_anchor])

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
    parser.add_argument("-trackpy_percentile", type=int, default=50)
    parser.add_argument("-known_anchor", type=str, default="c01 anchor")
    parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-chunk_size", type=int, default=4000)

    args = parser.parse_args()

    main(args)
