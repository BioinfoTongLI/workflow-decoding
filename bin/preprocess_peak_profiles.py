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
import logging
try:
    import cupy as xp
    logging.info(cp.__version__)
    cp.cuda.Device(0).use()
except ImportError:
    import numpy as xp
import pandas as pd

def main(stem, profiles, spot_loc, cleanup:bool=False):
    spot_profile = xp.load(profiles) if xp == cupy else np.load(profiles)
    spot_loc = pd.read_csv(spot_loc)
    if cleanup:
        xp_profile = xp.array(np.nan_to_num(spot_profile))
        # remove any spot that has  zero intensity in any channel/cycle
        try:
            mask = (xp_profile > 0).any((1, 2)).get()
        except AttributeError:
            mask = (xp_profile > 0).any((1, 2))
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
