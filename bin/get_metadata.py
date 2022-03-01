#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Read metadata for decoding
"""
import argparse
from reading_data_functions import read_taglist_and_channel_info
from pandas import read_csv
import numpy as np
import pickle


def main(args):
    # read channel_info.csv and taglist.csv
    barcodes_01, K, R, C, gene_names, channels_info = read_taglist_and_channel_info(
        args.auxillary_file_dir,
        taglist_name=args.taglist_name,
        channel_info_name=args.channel_info_name,
    )

    np.save("barcodes_01.npy", barcodes_01)
    np.save("gene_names.npy", gene_names)
    channels_info["K"] = K
    channels_info["R"] = R
    channels_info["C"] = C
    with open("channel_info.pickle", "wb") as fp:
        pickle.dump(channels_info, fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-auxillary_file_dir", type=str, required=True)
    parser.add_argument("-taglist_name", type=str, required=True)
    parser.add_argument("-channel_info_name", type=str, required=True)

    args = parser.parse_args()

    main(args)
