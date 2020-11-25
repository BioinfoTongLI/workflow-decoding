#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Load tiles and detects spots inside them
"""
import argparse
from spot_detection_functions import load_tiles_to_extract_spots
from pandas import read_csv
import numpy as np
import pickle


def main(args):

    # input parameters for loading the tiles and for spot detection
    tile_names = read_csv(args.tile_name)
    with open(args.channel_info, 'rb') as fp:
        channels_info = pickle.load(fp)
    C = channels_info["C"]
    R = channels_info["R"]

    tiles_to_load={'y_start':7, 'y_end':8, 'x_start':6, 'x_end':8} # select which tiles to load, including indices at the end
    tiles_info={'tile_size':args.tile_size, 'y_max':16, 'x_max':23, 'y_max_size':1000, 'x_max_size':1000, 'filename_prefix':args.filename_prefix}
    spots_params={'trackpy_diam_detect':5, 'trackpy_search_range':3, 'spot_diam_tophat':5, 'trackpy_prc':args.trackpy_percentile} # parameters for spot detection

    # load tile by tile and extract spots in each using trackpy
    spots, spots_loc, _ = load_tiles_to_extract_spots(args.tifs_dir, channels_info, C, R, tile_names, tiles_info, tiles_to_load, spots_params, args.tile_name_pattern, anchor_available=bool(args.anchor_available), coding_ch_starts_from=args.coding_ch_starts_from)

    # loading / saving already extracted spots
    # with open('%s_spots_trackpy_params.pickle' %args.stem, 'wb') as fp:
       # pickle.dump(spots_params, fp)
    np.save('%s_spots_trackpy_spot_profile.npy' %args.stem, spots)
    spots_loc.to_csv('%s_spots_trackpy_locations.csv' %args.stem, index=False)

    print('In total {} spots detected'.format(spots.shape[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-tifs_dir", type=str,
            required=True)
    parser.add_argument("-stem", type=str,
            required=True)
    parser.add_argument("-tile_name_pattern", type=str,
            required=True)
    parser.add_argument("-anchor_available", type=int,
            required=True)
    parser.add_argument("-tile_name", type=str,
            required=True)
    parser.add_argument("-tile_size", type=int,
            required=True)
    parser.add_argument("-filename_prefix", type=str,
            required=True)
    parser.add_argument("-trackpy_percentile", type=int,
            required=True)
    parser.add_argument("-channel_info", type=str,
            required=True)
    parser.add_argument("-coding_ch_starts_from", type=int,
            required=True)

    args = parser.parse_args()

    main(args)
