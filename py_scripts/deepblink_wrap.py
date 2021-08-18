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
import deepblink as db
import numpy as np
import skimage.io
import skimage.util
from pathlib import Path
from ome_zarr.reader import Reader
from ome_zarr.io import parse_url
import pandas as pd

import trackpy as tp


def deepblink_on_large(img, model, stepsize=1000):
    # Tile image into 1000x1000 patches with steps of 1000 to avoid multiple detections of the same spots
    windows = skimage.util.view_as_windows(
        img, window_shape=(1000, 1000), step=stepsize
    )
    preds = []
    for n_row, row in enumerate(windows):
        for n_col, col in enumerate(row):
            # Predict current window
            raw_pred = db.inference.predict(col, model).T
            # Add position back in large image
            pred = np.array(
                [raw_pred[0] + (n_row * stepsize), raw_pred[1] + (n_col * stepsize)]
            ).T
            preds.append(pred)
    coords = np.concatenate(preds)
    return coords


def main(zarr_in, stem, tpy_search_range):
    zarr_in = Path(zarr_in)
    assert zarr_in.exists()
    reader = Reader(parse_url(zarr_in))
    anchor_image = list(reader())[0].data[0]

    # Load model of deepblink for peak-calling
    model = db.io.load_model("/vesicle.h5")

    coords = deepblink_on_large(
        anchor_image.squeeze().max(axis=0).compute(),
        model,
    )

    df = pd.DataFrame(coords, columns=["y", "x"])
    df = df.assign(
        x_int=lambda df: np.round(df.x).astype(np.uint32),
        y_int=lambda df: np.round(df.y).astype(np.uint32),
    )
    df.to_csv(f"{stem}_max_projected_peaks.tsv", sep="\t")
    # all_coords = []
    # for i in range(anchor_image.squeeze().shape[0]):
    # coords = deepblink_on_large(anchor_image.squeeze()[i])
    # df = pd.DataFrame(coords, columns=["y", "x"])
    # df["frame"] = i
    # all_coords.append(df)
    # combine = pd.concat(all_coords)

    # # Track peaks across cycles
    # tracks = tp.link(combine, tpy_search_range, memory=0, adaptive_step=0.98)
    # tracks.to_csv(f"{stem}_tracked_peaks.tsv", sep="\t")


if __name__ == "__main__":
    fire.Fire(main)
