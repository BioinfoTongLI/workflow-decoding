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
import pandas as pd
import re
import numpy as np


channel_map = {"Cy5":"G", "AF488":"C", "Cy3":"T", "AF750":"A"}
nucleotide_map = {"1":"A", "2":"G", "3":"C", "4":"T"}

def convert2AGCT(codelist):
    AGCT_codes = []
    for g in codelist:
        AGCT_codes.append("".join([nucleotide_map[i] for i in str(g)]))
    return AGCT_codes


def get_channel_info(col_list):
    n_chs = []
    n_cycs = []
    ch_info = {}
    for col in col_list:
        m = re.search("cycle(\d+)_channel(\d+)_(.*)", col)
        if m:
            cyc = int(m.group(1))
            ch = int(m.group(2))
            n_chs.append(ch)
            n_cycs.append(cyc)
            ch_name = m.group(3)
            if ch_name == "DAPI":
                ch_info[ch_name] = "nuclei"
            else:
                ch_info[ch_name] = nucleotide_map[str(ch)]
            print(cyc, ch, ch_name, col)
    n_ch_set = set(n_chs)
    assert len(n_ch_set) == max(n_ch_set)
    n_cyc_set = set(n_cycs)
    assert len(n_cyc_set) == max(n_cyc_set)
    ch_info["nCycles"] = max(n_cyc_set)
    ch_info["nChannel"] = max(n_ch_set)
    return ch_info

def main(csv_file):
    d = pd.read_csv(csv_file)
    code_sizes = [len(str(c)) for c in d.code]
    n_cycle = np.unique(code_sizes)
    assert len(n_cycle) == 1
    channel_info = {}
    channel_info["nCycles"] = n_cycle[0]

    channel_dict = {}
    for col in d.columns:
        if col.startswith("cycle"):
            m = re.search("cycle(\d+)_channel(\d+)_(.*)", col)
            channel_dict[(m.group(1), m.group(2))] = m.group(3)
    print(channel_dict)
    nucleotide_codes = []
    for gene_ind in d.index:
        gene = d.loc[gene_ind]
        str_l = []
        for i, ind in enumerate(str(gene.code)):
            ch_name = channel_dict[(str(i+1), ind)]
            col_name = f"cycle{i+1}_channel{ind}_{ch_name}"
            assert gene[col_name] == 1
            str_l.append(ch_name)
        nucleotids = "".join([channel_map[s] for s in str_l])
        # print(nucleotids)
        nucleotide_codes.append(nucleotids)
    d["nucleotide_codes"] = nucleotide_codes
    channel_indexes = [int(k[1]) for k in channel_dict.keys()]
    assert len(np.unique(channel_indexes)) == np.max(channel_indexes)
    channel_info["nChannel"] = np.max(channel_indexes)
    channel_info["DAPI"] = "nuclei"
    for ch in channel_map:
        print(ch)
        if ch == "AF750":
            channel_info["Atto425"] = channel_map[ch]
        else:
            channel_info[ch] = channel_map[ch]
    print(channel_info)

    df = pd.DataFrame(pd.Series(channel_info)).T
    df.to_csv("channel_info.csv", index=False)
    print(df)

    # ch_info = get_channel_info(d.columns)
    # print(ch_info)
    # df = pd.DataFrame(pd.Series(ch_info)).T
    # columns_order = ["nCycles", "nChannel", "DAPI"]
    # nucleotides = [col for col in df.columns.values if col not in columns_order]
    # ordered_cols = columns_order + nucleotides
    # df = df[ordered_cols]

    # df.to_csv("channel_info.csv", index=False)
    # d["Channel"] = convert2AGCT(d.code)

    taglist = d[["gene", "nucleotide_codes"]]
    taglist = taglist.rename(columns={"gene":"Gene", "nucleotide_codes":"Channel"})
    print(taglist)
    pd.DataFrame(taglist).to_csv("taglist.csv", index=False)


if __name__ == "__main__":
    fire.Fire(main)
