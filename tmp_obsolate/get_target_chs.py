#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Preprocess images before GMM-decoding
"""
import argparse
import tifffile
from dask_image import imread
from apeer_ometiff_library import omexmlClass
from skimage.restoration import denoise_wavelet
import numpy as np
import re
import zarr

# from distributed import Client
# client = Client(n_workers=10, threads_per_worker=1)


def main(args):
    img = imread.imread(args.ome_tif)

    with tifffile.TiffFile(args.ome_tif) as tif:
        omexml_string = tif.ome_metadata

    metadata = omexmlClass.OMEXML(omexml_string)

    pixels = metadata.image(0).Pixels
    dim_format = pixels.DimensionOrder

    ch_names = [
        pixels.Channel(i).get_Name().replace(" ", "_") for i in range(pixels.SizeC)
    ]

    anchor_ch_indexes = []
    dic_cycle_ch = {}
    for i, n in enumerate(ch_names):
        current_c_index = int(re.search("c0(\d)_*", n).group(1))
        if current_c_index not in dic_cycle_ch:
            dic_cycle_ch[current_c_index] = []

        if n.endswith(args.anchor_ch_suffix):
            anchor_ch_indexes.append(i)
            continue

        if current_c_index > args.first_decoding_cycle and "DAPI" not in n:
            print(n)

            dic_cycle_ch[current_c_index].append(i)
            if n.endswith(args.anchor_ch_suffix):
                anchor_ch_indexes.append(i)
    print(dic_cycle_ch)
    print(anchor_ch_indexes)
    zarr.save("anchors.zarr", np.array(img[anchor_ch_indexes]))
    print(img[anchor_ch_indexes].shape)
    for cyc_i in dic_cycle_ch:
        if cyc_i > args.first_decoding_cycle:
            zarr.save("cycle_%s.zarr" % cyc_i, np.array(img[dic_cycle_ch[cyc_i]]))
            print(img[dic_cycle_ch[cyc_i]].shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-ome_tif", type=str, required=True)
    parser.add_argument("-first_decoding_cycle", type=int, required=True)
    parser.add_argument("-anchor_ch_suffix", type=str, default="c01_Alexa_647")

    args = parser.parse_args()

    main(args)
