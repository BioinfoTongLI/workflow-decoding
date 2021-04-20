#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

include { Extract_ch; Projection } from workflow.projectDir + '/nf_module_image_preprocessing/preprocess.nf'

params.ome_tif = 'path/to/ome.tiff'
params.out_dir = "./test/"
params.known_anchor = "c01 Alexa 647"
params.trackpy_separation = 2
params.rna_spot_size = 5
params.trackpy_percentile = 90

params.decode = true
params.auxillary_file_dir = "/nfs/team283_imaging/NT_ISS/playground_Tong/KR0018/new_opt/gmm-input/"
params.taglist_name = "taglist.csv"
params.channel_info_name = "channel_info.csv"

// not used in this version
params.coding_ch_starts_from = 0
params.anchor_available = 1

/*ome_tif_ch = Channel.fromPath(params.ome_tif).*/
    /*into{ome_tif_for_anchor_peak_calling; ome_tif_for_peak_intensity_extraction}*/

process Get_meatdata {
    echo true
    cache "lenient"
    /*storeDir params.out_dir + "decoding_metadata"*/
    publishDir params.out_dir + "decoding_metadata", mode:"copy"
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro"

    when:
    params.decode

    input:
    path gmm_input_dir
    val taglist_name
    val channel_info_name

    output:
    path "barcodes_01.npy", emit: barcodes
    path "gene_names.npy", emit: gene_names
    path "channel_info.pickle", emit: channel_infos

    script:
    """
    python /gmm_decoding/get_metadata.py -auxillary_file_dir ${gmm_input_dir}/  -taglist_name ${taglist_name} -channel_info_name ${channel_info_name}
    """
}

process Preprocess_anchor_image {
    echo true
    /*storeDir params.out_dir + "anchor_peaks"*/
    publishDir params.out_dir + "anchor_peaks", mode:"copy"
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro"

    input:
    file ome_tif from ome_tif_for_anchor_peak_calling

    output:
    tuple val(stem), file("${stem}_anchor.zarr") into processed_anchor_img_zarr

    script:
    stem = file(ome_tif).baseName
    """
    python /gmm_decoding/preprocess_anchor_chs.py -ome_tif ${ome_tif} -known_anchor "${params.known_anchor}" -stem ${stem}
    """
}

process Call_peaks_in_anchor {
    echo true
    /*storeDir params.out_dir + "anchor_peaks"*/
    publishDir params.out_dir + "anchor_peaks", mode:"copy"
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro"

    input:
    tuple stem, file(anchor_zarr) from processed_anchor_img_zarr

    output:
    tuple stem, file("${stem}_peaks.tsv") into processed_anchor

    script:
    """
    python /gmm_decoding/call_peaks_in_anchor.py -anchor_zarr ${anchor_zarr} -spot_diameter ${params.rna_spot_size} -trackpy_percentile ${params.trackpy_percentile} -trackpy_search_range 9 -peak_separation ${params.trackpy_separation} -stem ${stem}
    """
}

process Extract_coding_chs_to_zarr {
    echo true
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro"
    /*storeDir params.out_dir + "coding_chs"*/
    publishDir params.out_dir + "coding_chs", mode:"copy"

    when:
    params.decode

    input:
    file ch_info from channels_info
    file ome_tif from ome_tif_for_peak_intensity_extraction

    output:
    tuple val(stem),  file("${stem}.zarr") into decoding_chs

    script:
    stem = ome_tif.baseName
    """
    python /gmm_decoding/extract_coding_chs.py -ome_tif ${ome_tif} -ch_info ${ch_info} -out ${stem}.zarr
    """
}


process Preprocess_coding_chs {
    echo true
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro --gpus all"
    /*storeDir params.out_dir + "coding_chs_processed"*/
    publishDir params.out_dir + "coding_chs_processed", mode:"copy"

    when:
    params.decode

    input:
    tuple val(stem), file(raw_coding_chs) from decoding_chs

    output:
    tuple val(stem), file("${stem}_processed.zarr") into processed_chs

    script:
    """
    python /gmm_decoding/preprocess_coding_chs.py -zarr ${raw_coding_chs} -out ${stem}_processed.zarr -spot_diameter ${params.rna_spot_size}
    """
}


process Extract_peak_intensities_from_preprocessed_arrays {
    echo true
    /*storeDir params.out_dir + "peak_intensities"*/
    publishDir params.out_dir + "peak_intensities", mode:"copy"
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro --gpus all"

    when:
    params.decode

    input:
    tuple val(stem), file(anchor_peaks), file(zarr) from processed_anchor.combine(processed_chs, by:0)
    /*tuple val(stem), file(zarr) from */

    output:
    tuple val(stem), file("extracted_peak_intensities.npy"), file("spot_locs.csv") into paeks_for_decoding

    script:
    """
    python /gmm_decoding/extract_peak_intensities.py -zarr ${zarr} -anchors ${anchor_peaks} -out ${stem}
    """
}


process Decode_peaks {
    echo true
    /*storeDir params.out_dir + "decoded"*/
    publishDir params.out_dir + "decoded", mode:"copy"
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro --gpus all"

    when:
    params.decode

    input:
    tuple val(stem), file(spot_profile), file(spot_loc) from paeks_for_decoding
    file barcodes_f from barcodes_01
    file gene_names_f from gene_names
    file channel_info_f from channel_info_for_decode

    output:
    file "${stem}_decoded_df.tsv" into decoded_df
    file "${stem}_decode_out_parameters.pickle" into decode_out_parameters

    script:
    """
    python /gmm_decoding/decode.py -spot_profile ${spot_profile} -spot_loc ${spot_loc} -barcodes_01 ${barcodes_f} -gene_names ${gene_names_f} -channels_info ${channel_info_f} -stem ${stem}
    """
}


process Do_Plots {
    echo true
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro --gpus all"
    publishDir params.out_dir + "plots", mode:"copy"

    when:
    params.decode

    input:
    file decoded_df_f from decoded_df
    file decode_out_parameters_f from decode_out_parameters
    file channel_info_f from channel_info_for_plot

    output:
    file "*.png"

    script:
    """
    python /gmm_decoding/do_plots.py -decoded_df $decoded_df_f -decode_out_params $decode_out_parameters_f -channels_info ${channel_info_f}
    """
}


workflow {
    Get_meatdata(params.auxillary_file_dir, params.taglist_name, params.channel_info_name)
    Get_meatdata.out[1].view()
    println Get_meatdata.out[0]
    /*Extract_ch()*/
}