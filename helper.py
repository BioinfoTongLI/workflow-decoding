#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Python helper functions for GMM decoding
"""
import fire
from pathlib import Path
from ome_zarr.reader import Reader
from ome_zarr.writer import write_image
from ome_zarr.scale import Scaler
from ome_zarr.io import parse_url

# from cucim.skimage.morphology import white_tophat, disk
from skimage.morphology import white_tophat, disk
import trackpy as tp
import numpy as np
import pathlib
import zarr
import pickle


def white_tophat_cp(chunk, **kwargs):
    import cupy as cp

    return white_tophat(cp.asarray(chunk), **kwargs)


class Helper(object):
    def __init__(self, zarr_in: str):
        zarr_in = Path(zarr_in)
        assert zarr_in.exists()
        reader = Reader(parse_url(zarr_in))
        self.raw_data = list(reader())[0].data[0]

    def enhance_spots(self, stem: str, diam: int, ch_info: str, anchor_ch_ind=None):
        # from cucim.skimage.exposure import equalize_adapthist
        from skimage.exposure import equalize_adapthist, equalize_hist
        import dask.array as da
        from dask.distributed import Client, LocalCluster
        import tifffile as tf

        client = Client(
            n_workers=15,
            # processes=False,
            memory_limit="300GB",
        )
        print(client)

        with open(ch_info, "rb") as fp:
            channels_info = pickle.load(fp)
        n_ch_per_cycle = len(channels_info["channel_names"])
        n_cycle = channels_info["R"]
        if n_ch_per_cycle * n_cycle < self.raw_data.shape[1]:
            assert anchor_ch_ind
            first_cycle = [False] * n_ch_per_cycle
            first_cycle[anchor_ch_ind] = True
        # coding_ch_indexs = channels_info["coding_chs"] * channels_info["R"]
        projected_chs = []
        projected_chs.append(self.raw_data[0, anchor_ch_ind, 0])
        for i in range(channels_info["R"]):
            cyc_ch_indexes = (i + 1) * n_ch_per_cycle + np.where(
                channels_info["coding_chs"]
            )[0]
            projected_chs.append(np.max(self.raw_data[0, cyc_ch_indexes, 0], axis=0))
        chs_with_peaks = da.array(projected_chs)

        # This was intended to perform a normalization before peak enhancement,
        # however was having issues with dask and memory management
        # all_ch_with_peaks = all_ch_with_peaks.rechunk((1, 5 * 10**3, 5 * 10**3))
        # normed_chs_with_peaks = all_ch_with_peaks.map_blocks(
        # # equalize_adapthist, dtype=np.float16
        # equalize_hist, dtype=np.float64, nbins=180
        # )
        # normed_chs_with_peaks = normed_chs_with_peaks.rechunk((1, 2 ** 10, 2 ** 10))
        print(chs_with_peaks)
        hat_enhenced = chs_with_peaks.map_overlap(
            white_tophat,
            selem=np.expand_dims(disk(diam), 0),
            depth=(0, 100, 100),
            dtype=np.float16,
        )
        hat_enhenced = hat_enhenced.compute()
        store = parse_url(pathlib.Path(f"{stem}_spot_enhanced"), mode="w").store
        group = zarr.group(store=store).create_group("0")

        write_image(image=hat_enhenced, group=group, chunks=(2 ** 10, 2 ** 10))

        # tf.imwrite(f"{stem}_spot_enhanced.tif", hat_enhenced[0].squeeze(), imagej=True, metadata={'axes': 'YX'})

        client.close()
        return self

    def call_peaks(
        self,
        diam: int,
        stem: str,
        tp_percentile: int,
        peak_separation: int,
        tpy_search_range: int,
    ):
        df = tp.batch(
            self.raw_data.squeeze().compute(),
            diam,
            separation=peak_separation,
            percentile=tp_percentile,
            minmass=50,
            engine="numba",
        )
        # df.to_csv(f"{stem}_detected_peaks.tsv", sep="\t")
        t = tp.link(df, tpy_search_range, memory=0)

        tracks = tp.filtering.filter_stubs(t, 5)
        tracks = tracks.assign(
            x_int=lambda df: np.round(df.x).astype(np.uint64),
            y_int=lambda df: np.round(df.y).astype(np.uint64),
        )
        tracks.to_csv(f"{stem}_tracked_peaks.tsv", sep="\t")


if __name__ == "__main__":
    fire.Fire(Helper)
