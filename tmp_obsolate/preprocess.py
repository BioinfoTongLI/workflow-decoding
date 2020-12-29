#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Preprocess all channels for decoding
"""
import argparse
import numpy as np
import zarr
from skimage.restoration import denoise_wavelet
from skimage.morphology import white_tophat, disk
import dask.array as da
from pathlib import Path
import itk
import pickle
import cv2 as cv


def downsample(image, rescale_factor):
    small = cv.resize(
        image,
        (int(image.shape[1] / rescale_factor), int(image.shape[0] / rescale_factor)),
    )

    image_ikt = itk.image_view_from_array(small)
    image_ikt.SetSpacing([rescale_factor, rescale_factor])
    return image_ikt


def main(args):
    # xmin = 0
    # xmax = 9000
    # ymin = 0
    # ymax = 9000
    anchors = zarr.load(args.anchor_zarr).squeeze()
    # anchors = zarr.load("anchors.zarr").astype(np.float32).squeeze()[ymin:ymax, xmin:xmax]
    print(anchors.shape)

    img = zarr.load(args.decode_zarr)
    zproj = da.nanmax(da.from_array(img), axis=0).compute()
    print(type(zproj))
    print(zproj.shape)
    rescale_factor = 10
    anchor_small = downsample(anchors.astype(np.float32), rescale_factor)
    zproj_small = downsample(zproj.astype(np.float32), rescale_factor)
    print(anchor_small.shape)
    print(zproj_small.shape)

    parameter_object = itk.ParameterObject.New()
    parameter_map_rigid = parameter_object.GetDefaultParameterMap("bspline")
    parameter_object.AddParameterMap(parameter_map_rigid)

    # Call registration function
    result_image_small, result_transform_parameters = itk.elastix_registration_method(
        anchor_small, zproj_small, parameter_object=parameter_object
    )

    print(result_image_small.shape)
    stem = Path(args.decode_zarr).stem
    zarr.save("%s_registered.zarr" % stem, np.array(result_image_small))
    # with open('%s_registered.pickle' %stem, 'wb') as fp:
    # pickle.dump(params, fp)

    # img = da.array(img)
    # img = da.transpose(img, (1, 2, 0))
    # img = img.rechunk((2000, 2000, 1)).astype(np.float32)
    # print(img.shape)
    # print(img)

    # hat_enhenced = da.map_blocks(white_tophat, img,
    # selem=np.expand_dims(disk(args.whitehar_disk_diam), -1),
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

    parser.add_argument("-anchor_zarr", type=str, required=True)
    parser.add_argument("-decode_zarr", type=str, required=True)

    parser.add_argument("-overlaps", type=int, default=50)
    parser.add_argument("-whitehar_disk_diam", type=int, default=5)

    args = parser.parse_args()

    main(args)
