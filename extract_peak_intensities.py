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
import argparse
import numpy as np
from skimage.restoration import denoise_wavelet
from skimage.morphology import white_tophat, disk
import dask.array as da
import dask.dataframe as dd
from dask import delayed
from dask_image import imread
import tifffile as tf
from apeer_ometiff_library import omexmlClass
import pysnooper
import pickle


def make_in_range(l, upper):
    l = np.where(l < 0, 0, l)
    l = np.where(l > upper, upper, l)
    return l


# @delayed
def get_intensities(img, anchor_coords, r):
    img = np.array(img)
    shape = img.shape
    intensities = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            ys = anchor_coords[:, 0] + dy
            xs = anchor_coords[:, 1] + dx
            ys = make_in_range(ys, shape[0])
            xs = make_in_range(xs, shape[1])
            intensities.append(img[ys, xs, :, :])
            # print(img[ys, xs, :, :])
    return da.array(intensities)


@pysnooper.snoop()
def main(args):
    spots = dd.read_csv(args.anchors.name, sep="\t")[["y", "x"]].astype(np.uint32)
    anchor_coords = np.array(spots)
    imgs = da.from_zarr(args.zarr, "normalized_coding_chs")

    intensites = get_intensities(imgs, anchor_coords, args.r)

    np.save("extracted_peak_intensities.npy", da.max(intensites, axis=0).compute())
    spots = spots.compute()
    spots.columns = map(str.capitalize, spots.columns)
    spots["Tile"] = 0
    spots.to_csv("spot_locs.csv", index=False)

    # rescale_factor = 10
    # anchor_small = downsample(anchors.astype(np.float32), rescale_factor)
    # zproj_small = downsample(zproj.astype(np.float32), rescale_factor)
    # print(anchor_small.shape)
    # print(zproj_small.shape)

    # parameter_object = itk.ParameterObject.New()
    # parameter_map_rigid = parameter_object.GetDefaultParameterMap("bspline")
    # parameter_object.AddParameterMap(parameter_map_rigid)

    # # Call registration function
    # result_image_small, result_transform_parameters = itk.elastix_registration_method(
    # anchor_small, zproj_small, parameter_object=parameter_object
    # )

    # print(result_image_small.shape)
    # stem = Path(args.decode_zarr).stem
    # zarr.save("%s_registered.zarr" % stem, np.array(result_image_small))

    # with open('%s_registered.pickle' %stem, 'wb') as fp:
    # pickle.dump(params, fp)

    # img = da.array(img)
    # img = da.transpose(img, (1, 2, 0))
    # img = img.rechunk((2000, 2000, 1)).astype(np.float32)
    # print(img.shape)
    # print(img)

    # hat_enhenced = da.map_blocks(white_tophat, img,
    # selem=np.expand_dims(disk(args.whitehat_disk_diam), -1),
    # dtype=np.float32)
    # # print(hat_enhenced)
    # # denoised = da.map_overlap(
    # # denoise_wavelet, hat_enhenced,
    # # depth=(args.overlaps, args.overlaps, 0), trim=True,
    # # boundary="none",
    # # method='BayesShrink',
    # # mode='soft',
    # # multichannel=True).astype(np.uint16)
    # # print(denoised.compute())
    # hat_enhenced = da.transpose(hat_enhenced, (2, 0, 1))
    # zarr.save("%s_preprocessed.zarr" %Path(zarr_img).stem,
    # hat_enhenced.compute())

    # if len(anchor_ch_indexes) <= 1:
    # # create fake anchor channel
    # max_projs = []
    # for i in dic_cycle_ch:
    # if i > start_decode_from:
    # max_projs.append(np.max(img[dic_cycle_ch[i]], axis=0))
    # print(max_projs)
    # denoised_max_projs =[proj.map_blocks(denoise_wavelet, method='BayesShrink', mode='soft') for proj in max_projs]
    # print(denoised_max_projs)
    # im_bayes = [job.compute() for job in denoised_max_projs]
    # print(im_bayes)
    # else:
    # print(img[anchor_ch_indexes].shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-zarr", type=str, required=True)
    parser.add_argument("-anchors", type=argparse.FileType("r"), required=True)
    parser.add_argument("-out", type=str, required=True)
    parser.add_argument("-r", type=int, default=2)

    args = parser.parse_args()

    main(args)
