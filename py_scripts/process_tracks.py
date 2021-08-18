#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import cudf

# import dask_cudf
# import dask.dataframe as dd
import numpy as np


def complete_missing_cycles(df, n_cycle=6):
    part_id = df.iloc[0].loc["particle"]
    print(part_id)
    df = df.reset_index(drop=True).set_index("frame")
    for i in range(n_cycle):
        if i not in df.index:
            df = df.append(
                cudf.DataFrame(
                    [-1, None, None, part_id], index=df.columns, columns=[i]
                ).transpose()
            )
    df.sort_index(inplace=True)
    df.fillna(method="ffill", axis=0, inplace=True)
    return df


def main(tsv, stem):
    tracks = cudf.read_csv(tsv, sep="\t")
    # tracks = dd.read_csv(tsv, sep="\t")
    tracks = tracks.rename({"Unnamed: 0": "ID"}, axis=1)
    print(tracks.particle.unique(), tracks)
    # full_df = tracks.groupby("particle").apply(
    # complete_missing_cycles,
    # ).compute()
    # full_df = full_df.assign(
    # x_int=lambda df: np.round(df.x).astype(np.uint64),
    # y_int=lambda df: np.round(df.y).astype(np.uint64),
    # )
    # full_df.to_csv(f'{stem}_tracked_peaks.tsv', sep="\t")


if __name__ == "__main__":
    fire.Fire(main)
