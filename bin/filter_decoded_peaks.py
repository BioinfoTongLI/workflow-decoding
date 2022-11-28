#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the BSD-3 license.

"""
Filter GMM decoded spots
"""
import pandas as pd
import fire


def main(stem, transcripts, prob_threshold):
    decoded_spots_df = pd.read_csv(transcripts, sep="\t")
    selected_peaks = decoded_spots_df[decoded_spots_df.Probability > prob_threshold]
    selected_peaks.to_csv(
        f"{stem}_decoded_df_prob_thresholded_{prob_threshold}.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    fire.Fire(main)
