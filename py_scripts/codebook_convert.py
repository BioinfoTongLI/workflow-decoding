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
                ch_info[ch_name] = ch
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
    print(d.columns)
    ch_info = get_channel_info(d.columns)
    print(ch_info)
    pd.Series(ch_info).T.to_csv("channel_info.csv")
    d["Channel"] = convert2AGCT(d.code)
    taglist = d[["gene", "Channel"]]
    taglist = taglist.rename(columns={"gene":"Gene"})
    print(taglist)
    pd.DataFrame(taglist).to_csv("taglist.csv", index=False)


if __name__ == "__main__":
    fire.Fire(main)
