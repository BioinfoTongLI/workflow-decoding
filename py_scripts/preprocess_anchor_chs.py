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
import trackpy as tp
from aicsimageio import AICSImage


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
            sigma=spot_diameter * 1.5,
            rescale_sigma=True,
            multichannel=False,
        )
        * 10 ** 4
    ).astype(np.uint16)

    return denoised.compute()


def main(args):
    # Retrieve anchor channel
    imgs = AICSImage(args.ome_tif)

    print(imgs.dims, imgs.shape)
    print(imgs.get_image_dask_data("YX", T=0, Z=0, S=0, C=args.known_anchor_index))
    # processed_anchor = get_processed_anchor_ch(
    # imgs[args.known_anchor_index], args.chunk_size, args.spot_diameter
    # )

    # # Save
    # da.array(processed_anchor).to_zarr("%s_anchor.zarr" % args.stem, "processed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-ome_tif", type=str, required=True)
    parser.add_argument("-stem", type=str, required=True)
    parser.add_argument("-spot_diameter", type=int, default=5)
    parser.add_argument("-known_anchor_index", type=int, default=4)
    # parser.add_argument("-whitehat_disk_diam", type=int, default=5)
    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-chunk_size", type=int, default=4000)

    args = parser.parse_args()

    main(args)
