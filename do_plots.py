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


def plot_loss(out):
    print('The initial loss is {} and the final loss is {}'.format(1/out['class_probs'].shape[0]*out['params']['losses'][0],1/out['class_probs'].shape[0]*out['params']['losses'][len(out['params']['losses'])-1]))
    plt.figure(num=None, figsize=(4, 2), dpi=100, facecolor='w', edgecolor='k')
    plt.plot(np.arange(10,len(out['params']['losses'])),(1/out['class_probs'].shape[0]*np.asarray(out['params']['losses'][10:])))
    plt.title('Loss over iterations')
    plt.show()


def plot_mean_cov_of_classes(out, R, C):
    # plot estimated parameters for class means and covariance
    plt.figure(num=None, figsize=(16, 3), dpi=100, facecolor='w', edgecolor='k')
    activation = (out['params']['codes_tr_v_star']+out['params']['codes_tr_consts_v_star'])[0,:].numpy() #corresponding to the channel activation (code=1)
    no_activation = out['params']['codes_tr_consts_v_star'][0,:].numpy() # (code=0)
    channel_activation=np.stack((no_activation,activation))
    plt.subplot(1, 2, 1)
    plt.scatter(np.arange(1,1+R*C),activation,c='green')
    plt.scatter(np.arange(1,1+R*C),no_activation,c='orange')
    plt.legend(('channel active','not active'),bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.vlines(np.arange(0.5,R*C+.8,C), out['params']['codes_tr_consts_v_star'].min(), (out['params']['codes_tr_v_star']+out['params']['codes_tr_consts_v_star']).max(), linestyles='dashed')
    plt.title('Parameters of the barcode transformation as activation / no activation')
    plt.subplot(1, 2, 2)
    plt.imshow(out['params']['sigma_star'])
    plt.yticks(np.arange(3,R*C,4),np.arange(4, R*C+1, 4))
    plt.xticks(np.arange(3,R*C,4),np.arange(4, R*C+1, 4))
    plt.colorbar()
    plt.title('Estimated class covariance')
    #plt.savefig(os.getcwd() + '/out_imgs/' + dataset_name +'_params.png')
    plt.show()


def plot_hist_after_thresholding(decoded_df):
    # plot a histogram of class assignments when class probabilities are thresholded by thr
    thr=0.7
    df = pd.concat([decoded_df.Name[decoded_df.Probability > thr].replace('perforin','Perforin').value_counts(), decoded_df.Name[decoded_df.Probability <= thr].replace(np.unique(decoded_df.Name),'thr').value_counts()]).sort_index(axis=0)
    plt.figure(num=None, figsize=(14,4), dpi=100, facecolor='w', edgecolor='k')
    ax = df.plot(kind='bar',width=0.7,rot=70,logy=True,fontsize=6,figsize=(14,4))
    ax.set_facecolor('w')
    num_decoded_barcodes = sum((decoded_df.Name!='background')&(decoded_df.Name!='infeasible')&(decoded_df.Name!='NaN')&(decoded_df.Probability>thr))
    #ax.legend(["gmm: {} decoded spots with prob > {}".format(num_decoded_barcodes,thr)],fontsize=8)
    for p in ax.patches:
        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005),size=6)
    plt.title('Histogram of decoded barcodes afther thresholding with {}: \n in total {} spots detected while {} spots decoded ({}%)'.format(thr,decoded_df.shape[0], num_decoded_barcodes , np.round(100*num_decoded_barcodes/ decoded_df.shape[0], 2 )),fontsize=10)
    #plt.savefig(os.getcwd() + '/out_imgs/' + dataset_name +'_histogram.png')
    plt.show()

    print('Histogram of decoded barcodes afther thresholding with {}: \n in total {} spots detected while {} spots decoded ({}%)'.format(thr,decoded_df.shape[0], num_decoded_barcodes , np.round(100*num_decoded_barcodes/ decoded_df.shape[0], 2 )))

    print('Class names: {}'.format(np.unique(decoded_df.Name)))


def main(args):
    decoded_df = pd.read_csv(args.decoded_df, sep="\t")
    decoded_df["Name"] = decoded_df["Name"].astype(str)
    print(decoded_df.Probability)
    with open(args.decode_out_params, 'rb') as fp:
        decode_out_params = pickle.load(fp)
    # print(decode_out_params)
    with open(args.channels_info, 'rb') as fp:
        channels_info = pickle.load(fp)
    C = channels_info["C"]
    R = channels_info["R"]

    plot_loss(decode_out_params)
    plot_mean_cov_of_classes(decode_out_params, R, C)
    plot_hist_after_thresholding(decoded_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-decoded_df", type=str,
            required=True)
    parser.add_argument("-decode_out_params", type=str,
            required=True)
    parser.add_argument("-channels_info", type=str,
            required=True)

    args = parser.parse_args()

    main(args)
