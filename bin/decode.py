#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Take previously detected spots and GMM decode them
"""
import argparse
import numpy as np
import pandas as pd
from decoding_functions import decoding_function, decoding_output_to_dataframe
import pickle
import fire


def main(stem, spot_profile, spot_loc, barcodes_01, gene_names, channels_info):
    # load
    spot_profile = np.load(spot_profile)
    print(spot_profile)
    spot_loc = pd.read_csv(spot_loc)
    barcodes_01 = np.load(barcodes_01, allow_pickle=True)
    gene_names = np.load(gene_names, allow_pickle=True)
    with open(channels_info, "rb") as fp:
        channels_info = pickle.load(fp)

    # print(spot_profile.shape, spot_loc)
    # estimate GMM parameters and compute class probabilities
    out = decoding_function(spot_profile, barcodes_01, print_training_progress=True)
    # creating a data frame from the decoding output
    df_class_names = np.concatenate((gene_names, ["infeasible", "background", "nan"]))
    df_class_codes = np.concatenate(
        (channels_info["barcodes_AGCT"], ["NA", "0000", "NA"])
    )
    decoded_spots_df = decoding_output_to_dataframe(out, df_class_names, df_class_codes)
    decoded_df = pd.concat([decoded_spots_df, spot_loc], axis=1)

    decoded_df.to_csv(f"{stem}_decoded_df.tsv", sep="\t", index=False)

    with open(f"{stem}_decode_out_parameters.pickle", "wb") as fp:
        pickle.dump(out, fp, protocol=4)


if __name__ == "__main__":
    fire.Fire(main)
