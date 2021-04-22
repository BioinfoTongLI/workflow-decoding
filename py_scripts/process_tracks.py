#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Process tracks of anchor images
"""
import argparse
import trackpy as tp
# import dask.dataframe as dd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def main(args):
    tracks = pd.read_csv(args.tracks, sep="\t")
    print(len(tracks.particle.unique()))
    tracks = tp.filtering.filter_stubs(tracks, 2)
    print(len(tracks.particle.unique()))
    tracks = tracks.assign(
            x_int=lambda df: np.round(df.x).astype(np.uint64),
            y_int=lambda df: np.round(df.y).astype(np.uint64),
            )
    print(tracks)
    tracks.to_csv('%s_processed_peaks.tsv' %args.stem, sep="\t")
    # tp.plot_traj(tracks)
    # plt.savefig("trajs.png")
    # plt.close()
    # d = tp.compute_drift(tracks)
    # d.plot()
    # plt.savefig("drift.png")
    # plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-tracks", type=str, required=True)
    parser.add_argument("-stem", type=str, required=True)

    args = parser.parse_args()

    main(args)
