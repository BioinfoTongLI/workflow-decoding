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


def map_to_index(peak_row):
    indexes = []
    NAME_INDEX_MAP = {"A": "4", "T": "3", "G": "2", "C": "1"}
    for n in peak_row.Code:
        if (
            not isinstance(n, float)
            and n != "infeasible"
            and n != "background"
            and n != "0000"
            and n.lower() != "nan"
            and n != "NA"
        ):
            indexes.append("".join([NAME_INDEX_MAP[i] for i in n]))
        else:
            indexes.append("")
    return indexes


def get_column_names(n_ch, n_cyc):
    R_C_index = []
    for j in range(n_ch):
        for i in range(n_cyc):
            R_C_index.append(f"R{i}_C{j}")
    return R_C_index


def decode(
    stem,
    spot_profile,
    spot_loc,
    barcodes_01,
    gene_names,
    channels_info,
    chunk_size=3 * 10**6,
    n_cycle=6,
    n_ch=4,
):
    R_C_index = get_column_names(n_ch, n_cycle)
    # load
    spot_profile = np.load(spot_profile, allow_pickle=True)
    profile_df = pd.DataFrame(
        spot_profile.reshape(spot_profile.shape[0], -1), columns=R_C_index
    )

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
        decoded_spots_df = pd.concat(
            [
                decoding_output_to_dataframe(decoded, df_class_names, df_class_codes)
                for decoded in decoded_list
            ]
        )
    else:
        # estimate GMM parameters and compute class probabilities
        out = decoding_function(
            spot_profile, barcodes_01, print_training_progress=False
        )
        # creating a data frame from the decoding output
        decoded_spots_df = decoding_output_to_dataframe(
            out, df_class_names, df_class_codes
        )
        with open(f"{stem}_decode_out_parameters.pickle", "wb") as fp:
            pickle.dump(out, fp, protocol=4)

    assert decoded_spots_df.shape[0] == spot_loc.shape[0]
    decoded_spots_df.loc[:, "y_int"] = spot_loc.y_int.values
    decoded_spots_df.loc[:, "x_int"] = spot_loc.x_int.values
    decoded_spots_df = decoded_spots_df.assign(index_code=map_to_index, axis=1)
    decoded_spots_df = pd.concat([decoded_spots_df, profile_df], axis=1)
    # assert decoded_spots_df.isnull().values.any() # shouldn't have any nan in the df

    decoded_spots_df.to_csv(f"{stem}_decoded_df.tsv", sep="\t", index=False)


if __name__ == "__main__":
    import fire
    fire.Fire(decode)
