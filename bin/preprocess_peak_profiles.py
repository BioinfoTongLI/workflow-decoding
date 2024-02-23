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
try:
    import cupy as xp
except ImportError:
    import numpy as xp
import pandas as pd

def main(stem, profiles, spot_loc, cleanup:bool=False):
    spot_profile = np.load(profiles)
    spot_loc = pd.read_csv(spot_loc)
    if cleanup:
        xp_profile = xp.array(np.nan_to_num(spot_profile))
        mask = (xp_profile > 0).any((1, 2)).get()  # remove any spot that has  zero intensity in any channel/cycle
        filtered_profiles = spot_profile[mask]
        filtered_spot = spot_loc[mask]
    else:
        filtered_profiles = np.nan_to_num(spot_profile)
        filtered_spot = spot_loc
    np.save(
        f"{stem}_filtered_peak_intensities.npy",
        filtered_profiles,
        allow_pickle=True,
    )
    filtered_spot.to_csv(f"{stem}_filtered_peak_locs.csv")


if __name__ == "__main__":
    fire.Fire(main)
