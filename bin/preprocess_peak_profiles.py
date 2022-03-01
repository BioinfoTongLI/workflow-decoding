#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import numpy as np
import cupy as cp
import pandas as pd


def main(profiles, spot_loc, stem):
    spot_profile = np.load(profiles)
    # print(spot_profile)
    spot_loc = pd.read_csv(spot_loc)
    cp_profile = cp.array(spot_profile)
    mask = [0.0 not in p.ravel() for p in cp_profile]
    filtered_profiles = spot_profile[mask]
    np.save(
        f"{stem}_filtaered_peak_intensities.npy",
        filtered_profiles,
        allow_pickle=True,
    )
    filtered_spot = spot_loc[mask]
    filtered_spot.to_csv(f"{stem}_filtered_peak_locs.csv")


if __name__ == "__main__":
    fire.Fire(main)
