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

try:
    from cucim.skimage.morphology import white_tophat, disk
    # from cucim.skimage.exposure import rescale_intensity
    import cupy as xp
    print("Using cucim")
except:
    from skimage.morphology import white_tophat, disk
    # from skimage.exposure import rescale_intensity
    # from skimage.restoration import denoise_wavelet
    import numpy as xp
    print("Using skimage")
import trackpy as tp
import numpy as np
import pathlib
import zarr
import pickle
import dask.array as da


def white_tophat_xp(chunk, **kwargs):
    try:
        return white_tophat(xp.array(chunk), **kwargs).get().astype(np.uint16)
    except:
        return white_tophat(xp.array(chunk), **kwargs).astype(np.uint16)


class Helper(object):
    def __init__(self, zarr_in: str):
        zarr_in = Path(zarr_in)
        assert zarr_in.exists()
        reader = Reader(parse_url(zarr_in))
        self.raw_data = list(reader())[0].data[0]

    def enhance_spots(
        self,
        stem: str,
        diam: int,
        ch_info: str,
        anchor_ch_ind: int = None,
        tophat: bool = False,
    ):
        # from cucim.skimage.exposure import equalize_adapthist
        # from skimage.exposure import equalize_adapthist, equalize_hist
        import tifffile as tf

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
        if tophat:
            hat_enhenced = chs_with_peaks.map_overlap(
                white_tophat,
                selem=np.expand_dims(disk(diam), 0),
                depth=(0, 20, 20),
                dtype=np.float16,
            )
            hat_enhenced = hat_enhenced.compute()
        else:
            hat_enhenced = chs_with_peaks.astype(np.float16).compute()
        store = parse_url(
            pathlib.Path(f"{stem}_spot_enhanced_diam_{diam}"), mode="w"
        ).store
        group = zarr.group(store=store).create_group("0")

        write_image(image=hat_enhenced, group=group, chunks=(2**10, 2**10))

        # tf.imwrite(f"{stem}_spot_enhanced.tif", hat_enhenced[0].squeeze(), imagej=True, metadata={'axes': 'YX'})

        return self

    def enhance_all(
            self, stem: str, diam: int, anchor_ch_index: int, whitehat: bool = False
    ):
        chs_with_peaks = self.raw_data[0, :, 0]
        print(chs_with_peaks)
        footprint = disk(diam // 2)

        store = parse_url(
            pathlib.Path(f"{stem}_spot_enhanced_diam_{diam}"), mode="w"
        ).store

        ref_img_to_match_hist = None
        for i in range(chs_with_peaks.shape[0]):
            ch = chs_with_peaks[i].rechunk({0: "auto", 1: "auto"})
            if whitehat:
                hat_enhenced = ch.map_overlap(
                    white_tophat_xp,
                    depth=(diam * 3, diam * 3),
                    footprint=footprint,
                    dtype=np.uint16,
                )
            else:
                hat_enhenced = ch

            # Wavelet denoise
            # denoised = (
                # hat_enhenced.map_overlap(
                    # denoise_wavelet,
                    # # depth=(args.overlaps, args.overlaps, 0, 0),
                    # depth=(diam * 3, diam * 3),
                    # method="BayesShrink",
                    # mode="soft",
                    # sigma=diam * 1.5,
                    # rescale_sigma=True,
                    # channel_axis=None,
                    # dtype=np.float16,
                # )
                # * 10**4
            # ).astype(np.uint16)

            # rescaled_ch = rescale_intensity(matched, in_range='image', out_range='dtype')
            group = zarr.group(store=store).create_group(f"0/{i}")
            write_image(image=hat_enhenced.compute(), group=group, axes="yx")

            del hat_enhenced
            try:
                xp._default_memory_pool.free_all_blocks()
            except:
                pass

    def call_peaks(
        self,
        diam: int,
        stem: str,
        tp_percentile: int,
        peak_separation: int,
        tp_search_range: int,
    ):
        print(self.raw_data)
        df = tp.locate(
            self.raw_data.compute().astype(np.uint16),
            diam,
            separation=peak_separation,
            percentile=tp_percentile,
            engine="numba",
            # minmass=50,
        )
        df["x_int"] = df.x.astype(np.uint32)
        df["y_int"] = df.y.astype(np.uint32)
        df.to_csv(
            f"{stem}_detected_peaks_diam_{diam}_percentile_{tp_percentile}_sep_{peak_separation}_search_range_{tp_search_range}.tsv",
            sep="\t",
            index=False,
        )
        # t = tp.link(df, tp_search_range, memory=0)

        # tracks = tp.filtering.filter_stubs(t, 5)
        # tracks = tracks.assign(
        # x_int=lambda df: np.round(df.x).astype(np.uint64),
        # y_int=lambda df: np.round(df.y).astype(np.uint64),
        # )
        # tracks.to_csv(f"{stem}_tracked_peaks.tsv", sep="\t")


if __name__ == "__main__":
    fire.Fire(Helper)
