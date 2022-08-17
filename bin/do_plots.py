#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Take decode output and do plots
"""
import argparse
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
import fire


def plot_loss(out):
    print(
        "The initial loss is {} and the final loss is {}".format(
            1 / out["class_probs"].shape[0] * out["params"]["losses"][0],
            1
            / out["class_probs"].shape[0]
            * out["params"]["losses"][len(out["params"]["losses"]) - 1],
        )
    )
    plt.figure(num=None, figsize=(4, 2), dpi=100, facecolor="w", edgecolor="k")
    plt.plot(
        np.arange(10, len(out["params"]["losses"])),
        (1 / out["class_probs"].shape[0] * np.asarray(out["params"]["losses"][10:])),
    )
    plt.title("Loss over iterations")
    # plt.show()
    plt.savefig("Loss.png")


def plot_mean_cov_of_classes(out, R, C):
    # plot estimated parameters for class means and covariance
    plt.figure(num=None, figsize=(16, 3), dpi=100, facecolor="w", edgecolor="k")
    activation = (
        out["params"]["codes_tr_v_star"] + out["params"]["codes_tr_consts_v_star"]
    )[
        0, :
    ].numpy()  # corresponding to the channel activation (code=1)
    no_activation = out["params"]["codes_tr_consts_v_star"][0, :].numpy()  # (code=0)
    channel_activation = np.stack((no_activation, activation))
    plt.subplot(1, 2, 1)
    plt.scatter(np.arange(1, 1 + R * C), activation, c="green")
    plt.scatter(np.arange(1, 1 + R * C), no_activation, c="orange")
    plt.legend(
        ("channel active", "not active"), bbox_to_anchor=(1.05, 1), loc="upper left"
    )
    plt.vlines(
        np.arange(0.5, R * C + 0.8, C),
        out["params"]["codes_tr_consts_v_star"].min(),
        (
            out["params"]["codes_tr_v_star"] + out["params"]["codes_tr_consts_v_star"]
        ).max(),
        linestyles="dashed",
    )
    plt.title("Parameters of the barcode transformation as activation / no activation")
    plt.subplot(1, 2, 2)
    plt.imshow(out["params"]["sigma_star"])
    plt.yticks(np.arange(3, R * C, 4), np.arange(4, R * C + 1, 4))
    plt.xticks(np.arange(3, R * C, 4), np.arange(4, R * C + 1, 4))
    plt.colorbar()
    plt.title("Estimated class covariance")
    plt.savefig("params.png")
    # plt.show()


def plot_hist_after_thresholding(decoded_df):
    # plot a histogram of class assignments when class probabilities are thresholded by thr
    thr = 0.7
    df = pd.concat(
        [
            decoded_df.Name[decoded_df.Probability > thr]
            .replace("perforin", "Perforin")
            .value_counts(),
            decoded_df.Name[decoded_df.Probability <= thr]
            .replace(np.unique(decoded_df.Name), "thr")
            .value_counts(),
        ]
    ).sort_index(axis=0)
    plt.figure(num=None, figsize=(14, 4), dpi=100, facecolor="w", edgecolor="k")
    ax = df.plot(kind="bar", width=0.7, rot=70, logy=True, fontsize=6, figsize=(14, 4))
    ax.set_facecolor("w")
    num_decoded_barcodes = sum(
        (decoded_df.Name != "background")
        & (decoded_df.Name != "infeasible")
        & (decoded_df.Name != "NaN")
        & (decoded_df.Probability > thr)
    )
    # ax.legend(["gmm: {} decoded spots with prob > {}".format(num_decoded_barcodes,thr)],fontsize=8)
    for p in ax.patches:
        ax.annotate(
            str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005), size=6
        )
    plt.title(
        "Histogram of decoded barcodes afther thresholding with {}: \n in total {} spots detected while {} spots decoded ({}%)".format(
            thr,
            decoded_df.shape[0],
            num_decoded_barcodes,
            np.round(100 * num_decoded_barcodes / decoded_df.shape[0], 2),
        ),
        fontsize=10,
    )
    plt.savefig("histogram.png")

    print(
        "Histogram of decoded barcodes afther thresholding with {}: \n in total {} spots detected while {} spots decoded ({}%)".format(
            thr,
            decoded_df.shape[0],
            num_decoded_barcodes,
            np.round(100 * num_decoded_barcodes / decoded_df.shape[0], 2),
        )
    )

    print("Class names: {}".format(np.unique(decoded_df.Name)))


# function creating a heatmap for plotting spatial patterns
def heatmap_pattern(decoded_df,name,grid=150, thr=0.7,plot_probs=True):
    if not 'Probability' in decoded_df.columns:
        if not 'Score' in decoded_df.columns:
            plot_probs = False
            x_coord = np.floor(decoded_df.X[(decoded_df.Name == name)].to_numpy(dtype=np.double) / grid).astype(np.int32)
            y_coord = np.floor(decoded_df.Y[(decoded_df.Name == name)].to_numpy(dtype=np.double) / grid).astype(np.int32)
        else:
            x_coord = np.floor(decoded_df.X[(decoded_df.Name == name) & (decoded_df.Score > thr)].to_numpy(dtype=np.double) / grid).astype(np.int32)
            y_coord = np.floor(decoded_df.Y[(decoded_df.Name == name) & (decoded_df.Score > thr)].to_numpy(dtype=np.double) / grid).astype(np.int32)
    else:
        x_coord = np.floor(decoded_df.X[(decoded_df.Name == name) & (decoded_df.Probability >thr)].to_numpy(dtype=np.double)/grid).astype(np.int32)
        y_coord = np.floor(decoded_df.Y[(decoded_df.Name == name) & (decoded_df.Probability >thr)].to_numpy(dtype=np.double)/grid).astype(np.int32)
    H = np.zeros((int(np.ceil(decoded_df.Y.to_numpy(dtype=np.double).max()/grid)),int(np.ceil(decoded_df.X.to_numpy(dtype=np.double).max()/grid))))
    if plot_probs:
        if 'Probability' in decoded_df.columns:
            prob = decoded_df.Probability[decoded_df.Name == name].to_numpy(dtype=np.double)
        elif 'Score' in decoded_df.columns:
            prob = decoded_df.Score[decoded_df.Name == name].to_numpy(dtype=np.double)
        prob[prob<thr]=0
        for i in range(len(x_coord)):
            H[y_coord[i],x_coord[i]] = H[y_coord[i],x_coord[i]] + prob[i]
    else:
        coords = np.concatenate((y_coord.reshape((len(x_coord),1)),x_coord.reshape((len(x_coord),1))), axis=1)
        coords_u ,coords_c = np.unique(coords ,axis=0, return_counts=True)
        H[coords_u[:,0],coords_u[:,1]]=coords_c
    return H


def main(decoded_df, decode_out_params):
    import os
    import sys
    import matplotlib.pyplot as plt
    sys.path.insert(1, os.path.join(sys.path[0], '../bin'))
    from do_plots import heatmap_pattern

    # with open(decode_out_params, "rb") as fp:
        # decode_out_params = pickle.load(fp)
    # print(decode_out_params)
    # with open(channels_info, "rb") as fp:
        # channels_info = pickle.load(fp)
    # C = channels_info["C"]
    # R = channels_info["R"]

    # plot_mean_cov_of_classes(decode_out_params, R, C)

    # plot_loss(decode_out_params)

    decoded_df = pd.read_csv(decoded_df, sep="\t")
    decoded_df["Name"] = decoded_df["Name"].astype(str)
    decoded_df["X"] = decoded_df.x_int
    decoded_df["Y"] = decoded_df.y_int


    plot_hist_after_thresholding(decoded_df)

    for n in decoded_df.Name.unique():
        hm = heatmap_pattern(decoded_df, n, grid=150)
        fig, ax = plt.subplots()
        pos = ax.imshow(hm)
        fig.colorbar(pos, ax=ax)
        plt.title(n)
        plt.savefig(f"./{n}.png", dpi=200)


if __name__ == "__main__":
    fire.Fire(main)
