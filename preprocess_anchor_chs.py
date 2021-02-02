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
import pysnooper


def get_processed_anchor_ch(anchor_img, chunk_size, spot_diameter):
    # anchor_img = da.transpose(anchor_img, (1, 2, 0))
    anchor_img = anchor_img.rechunk((1, chunk_size, chunk_size))

    # Tophat enhancement
    hat_enhenced = anchor_img.map_blocks(
        white_tophat,
        selem=np.expand_dims(disk(spot_diameter), 0),
        dtype=np.uint16,
    )

    # Wavelet denoise
    denoised = (
        hat_enhenced.map_overlap(
            denoise_wavelet,
            depth=(0, args.overlaps, args.overlaps),
            trim=True,
            boundary="nearest",
            method="BayesShrink",
            mode="soft",
            sigma=spot_diameter * 3,
            rescale_sigma=True,
            multichannel=False,
        )
        * 10 ** 4
    )

    return denoised.astype(np.uint16).compute()


@pysnooper.snoop()
def main(args):
    # Retrieve anchor channel
    imgs = imread.imread(args.ome_tif)
    with tf.TiffFile(args.ome_tif) as fh:
        metadata = omexmlClass.OMEXML(fh.ome_metadata)
    pixels = metadata.image(0).Pixels
    ch_names = np.array(pixels.get_channel_names())

    target_ch_masks = [args.known_anchor in ch for ch in ch_names]

    processed_anchor = get_processed_anchor_ch(
        imgs[target_ch_masks], args.chunk_size,
        args.spot_diameter
    )

    # Save
    da.array(processed_anchor).to_zarr("anchor_image.zarr", "processed_anchor")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-ome_tif", type=str, required=True)
    parser.add_argument("-spot_diameter", type=int, default=5)
    parser.add_argument("-known_anchor", type=str, default="c01 anchor")
    # parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-chunk_size", type=int, default=4000)

    args = parser.parse_args()

    main(args)
