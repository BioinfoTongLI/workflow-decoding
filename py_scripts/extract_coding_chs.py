#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Preprocess coding channels before peak-calling
"""
import argparse
import pysnooper
import pickle
import numpy as np
import tifffile as tf
from apeer_ometiff_library import omexmlClass

from dask_image import imread
import dask.array as da
from dask import delayed
from dask.distributed import Client
import zarr


def get_coding_ch_indexes(tif_path, ch_info, coding_cycle_starts_from):
    coding_ch_names = np.array(ch_info["channel_names"])[ch_info["coding_chs"]]
    with tf.TiffFile(args.ome_tif) as fh:
        metadata = omexmlClass.OMEXML(fh.ome_metadata)
    pixels = metadata.image(0).Pixels
    ch_names = np.array(pixels.get_channel_names())
    coding_ch_indexes = {}
    for i, ch in enumerate(ch_names):
        splits = ch.split(" ")
        dye_name = splits[-1]
        cycle = int(splits[0][-1])
        if (
            cycle >= coding_cycle_starts_from
            and dye_name in coding_ch_names
            and dye_name != "DAPI"
        ):
            if cycle not in coding_ch_indexes:
                coding_ch_indexes[cycle] = []
            coding_ch_indexes[cycle].append(i)
    return coding_ch_indexes


@pysnooper.snoop()
def main(args):
    # from distributed import Client
    # client = Client()
    with open(args.ch_info, "rb") as fp:
        channels_info = pickle.load(fp)
    C = channels_info["C"]
    R = channels_info["R"]

    coding_ch_indexes = get_coding_ch_indexes(
        args.ome_tif, channels_info, args.coding_cycle_starts_from
    )
    imgs = imread.imread(args.ome_tif)
    ch_base = np.array(channels_info["channel_base"])[channels_info["coding_chs"]]
    img_size = imgs.nbytes / 10 ** 9
    reshaped = []
    all_ch_imgs = da.asarray([imgs[coding_ch_indexes[ch]] for ch in coding_ch_indexes])
    reshaped_ch_imgs = da.transpose(all_ch_imgs, (2, 3, 1, 0))
    rechunked_ch_imgs = reshaped_ch_imgs.rechunk({0: "auto", 1: "auto", 2: 1, 3: 1})

    store = zarr.DirectoryStore(args.out)
    g = zarr.group(store=store)
    g["cycle"] = np.arange(1, 1 + R)
    g["base"] = ch_base
    rechunked_ch_imgs.to_zarr(args.out, "coding_ch_images")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-ome_tif", type=str, required=True)
    parser.add_argument("-ch_info", type=str, required=True)
    parser.add_argument("-out", type=str, required=True)

    parser.add_argument("-coding_cycle_starts_from", type=int, default=2)
    parser.add_argument("-chunk_size", type=int, default=2000)

    args = parser.parse_args()

    main(args)
