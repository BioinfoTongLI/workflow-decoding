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


def main(
    stem,
    spot_profile,
    spot_loc,
    barcodes_01,
    gene_names,
    channels_info,
    prob_threshold=0.6,
    chunk_size=3 * 10**6,
):
    # load
    spot_profile = np.load(spot_profile, allow_pickle=True)
    if spot_loc.endswith(".tsv"):
        spot_loc = pd.read_csv(spot_loc, index_col=0, sep="\t")
    else:
        spot_loc = pd.read_csv(spot_loc, index_col=0)
    print(spot_profile.shape, spot_loc)
    barcodes_01 = np.load(barcodes_01, allow_pickle=True)
    gene_names = np.load(gene_names, allow_pickle=True)
    with open(channels_info, "rb") as fp:
        channels_info = pickle.load(fp)
    df_class_codes = np.concatenate(
        (channels_info["barcodes_AGCT"], ["NA", "0000", "NA"])
    )
    df_class_names = np.concatenate((gene_names, ["infeasible", "background", "nan"]))
    print(df_class_names)

    n_spot = spot_profile.shape[0]
    n_chunk = int(np.round(n_spot / chunk_size)) + 1

    if n_spot >= chunk_size:
        # if there are too much spots
        decoded_list = [
            decoding_function(chunk, barcodes_01, print_training_progress=False)
            for chunk in np.array_split(spot_profile, n_chunk, axis=0)
        ]
        decoded_spots_df = pd.concat([decoding_output_to_dataframe(decoded, df_class_names, df_class_codes) for decoded in decoded_list])
    else:
        # estimate GMM parameters and compute class probabilities
        out = decoding_function(
            spot_profile, barcodes_01, print_training_progress=False
        )
        # creating a data frame from the decoding output
        decoded_spots_df = decoding_output_to_dataframe(out, df_class_names, df_class_codes)
        with open(f"{stem}_decode_out_parameters.pickle", "wb") as fp:
            pickle.dump(out, fp, protocol=4)

    assert decoded_spots_df.shape[0] == spot_loc.shape[0]
    decoded_spots_df.loc[:, "y_int"] = spot_loc.y_int.values
    decoded_spots_df.loc[:, "x_int"] = spot_loc.x_int.values
    # decoded_df = pd.concat([decoded_spots_df, spot_loc], axis=1)
    # assert decoded_spots_df.isnull().values.any() # shouldn't have any nan in the df

    decoded_spots_df.to_csv(f"{stem}_decoded_df.tsv", sep="\t", index=False)
    selected_peaks = decoded_spots_df[decoded_spots_df.Probability > prob_threshold]
    selected_peaks.to_csv(f"{stem}_decoded_df_prob_thresholded_{prob_threshold}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    fire.Fire(main)
