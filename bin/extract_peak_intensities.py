#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Use anchor coordinates to extract intensities
"""
import fire
import numpy as np
from skimage.restoration import denoise_wavelet
from skimage.morphology import white_tophat, disk
import dask.array as da
import dask.dataframe as dd
from dask_image import imread
import tifffile as tf

import trackpy as tp
import pickle
from ome_zarr.reader import Reader
from ome_zarr.io import parse_url
from pathlib import Path


def make_in_range(l, upper):
    l = np.where(l < 0, 0, l)
    l = np.where(l > upper, upper, l)
    return l


def get_intensities(img, peaks, r):
    # img = np.array(img)
    shape = img.shape
    intensities = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            ys = peaks[:, 0] + dy
            xs = peaks[:, 1] + dx
            ys = make_in_range(ys, shape[0])
            xs = make_in_range(xs, shape[1])
            intensities.append(img[ys, xs])
            # print(img[ys, xs])
    return da.array(intensities)


def main_multicycle(stem, raw_zarr, peaks, channel_info, coding_cyc_starts_from=1):
    peaks = dd.read_csv(peaks, sep="\t")[
        ["y_int", "x_int", "frame", "particle"]
    ].compute()
    peaks = peaks[peaks.frame >= coding_cyc_starts_from]
    peaks["frame"] = peaks.frame - coding_cyc_starts_from
    print(peaks)

    with open(channel_info, "rb") as fp:
        channel_info = pickle.load(fp)

    raw_zarr = Path(raw_zarr)
    assert raw_zarr.exists()
    reader = Reader(parse_url(raw_zarr))
    raw_data = list(reader())[0].data[0].squeeze()

    print(raw_data)
    # print(channel_info)
    R = channel_info["R"]
    peak_intensities = []
    for i in range(R):
        cyc_ch_indexes = (i + coding_cyc_starts_from) * len(
            channel_info["coding_chs"]
        ) + np.where(channel_info["coding_chs"])[0]
        print(cyc_ch_indexes)
        peaks_in_current_frame = peaks[peaks.frame == i]
        images_in_current_cycle = raw_data[cyc_ch_indexes, :, :].compute()
        curremt_peaks = images_in_current_cycle[
            :,
            np.array(peaks_in_current_frame.y_int),
            np.array(peaks_in_current_frame.x_int),
        ]
        curremt_peaks = np.transpose(curremt_peaks, (1, 0))
        print(curremt_peaks.shape)
        peak_intensities.append(curremt_peaks)
    print(np.array(peak_intensities).shape)
    # formatted_img = np.transpose(da.array(peak_intensities), (2, 3, 1, 0))
    print(formatted_img)

    # formatted_img = formatted_img.rechunk({0: -1, 1: -1, 2: 1, 3: 1})

    # print(formatted_img)

    # intensites = formatted_img.map_blocks(
    # get_intensities, peaks, 1, dtype=np.uint16
    # )

    # print(intensites.compute())

    # np.save(
    # f"{stem}_extracted_peak_intensities.npy",
    # da.max(intensites, axis=0).compute().astype(np.int16),
    # )
    # peaks = peaks.compute()
    # peaks.columns = map(str.capitalize, peaks.columns)
    # peaks["Tile"] = 0
    # peaks.to_csv(f"{stem}_peak_locs.csv", index=False)


def main(stem, raw_zarr, peaks, channel_info, coding_cyc_starts_from, peak_radius=1):
    peak_sur_coord = np.arange(-peak_radius, peak_radius + 1)
    xx, yy = np.meshgrid(peak_sur_coord, peak_sur_coord)
    dXYpair = zip(xx.reshape(-1), yy.reshape(-1))

    peaks = dd.read_csv(peaks, sep="\t")[["y_int", "x_int"]].compute()

    raw_zarr = Path(raw_zarr)
    assert raw_zarr.exists()
    reader = Reader(parse_url(raw_zarr))
    raw_data = list(reader())[0].data[0].squeeze()
    img_Y, img_X = raw_data.shape[-2:]
    print(img_Y, img_X)

    Xs, Ys = [], []
    for dx, dy in dXYpair:
        Ys.append(make_in_range(peaks.y_int + dy, img_Y))
        Xs.append(make_in_range(peaks.x_int + dx, img_X))

    with open(channel_info, "rb") as fp:
        channel_info = pickle.load(fp)
    # print(channel_info)

    R = channel_info["R"]
    peak_intensities = []
    for i in range(R):
        cyc_ch_indexes = (i + coding_cyc_starts_from) * len(
            channel_info["coding_chs"]
        ) + np.where(channel_info["coding_chs"])[0]
        print(cyc_ch_indexes)

        images_in_current_cycle = raw_data[cyc_ch_indexes, :, :].compute()
        surrounding_pixels = []
        for j in range(len(Xs)):
            surrounding_pixels.append(
                    images_in_current_cycle[:, Ys[j], Xs[j]]
                )
        surrounding_pixels = np.array(surrounding_pixels)
        max_intensities = np.max(surrounding_pixels, axis=0)
        print(max_intensities.shape)
        peak_intensities.append(np.transpose(max_intensities, (1, 0)))
    formatted_peak_profiles = np.transpose(
        da.array(peak_intensities), (1, 2, 0)
    ).astype(np.int32)
    # print(formatted_peak_profiles.shape)
    np.save(
        f"{stem}_extracted_peak_intensities.npy",
        formatted_peak_profiles,
        allow_pickle=True,
    )
    peaks.to_csv(f"{stem}_peak_locs.csv")


if __name__ == "__main__":
    fire.Fire(main)
