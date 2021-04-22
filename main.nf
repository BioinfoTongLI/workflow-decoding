#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

include { Extract_ch as extract_anchor_chs;
        Enhance_spots } from projectDir + '/nf_module_image_preprocessing/preprocess.nf'

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
    storeDir params.out_dir + "/decoding_metadata"
    /*publishDir params.out_dir + "/decoding_metadata", mode:"copy"*/

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
    python ${workflow.projectDir}/py_scripts/get_metadata.py -auxillary_file_dir ${gmm_input_dir}/  -taglist_name ${taglist_name} -channel_info_name ${channel_info_name}
    """
}

process Call_peaks_in_anchor {
    echo true
    storeDir params.out_dir + "/anchor_peaks"
    /*publishDir params.out_dir + "/anchor_peaks", mode:"copy"*/

    input:
    tuple val(stem), file(anchor_zarr)

    output:
    tuple val(stem), file("${stem}_peaks.tsv"), emit: peaks_from_anchor_ch

    script:
    """
    python ${workflow.projectDir}/py_scripts/call_peaks_in_anchor.py -anchor_zarr ${anchor_zarr} -spot_diameter ${params.rna_spot_size} -trackpy_percentile ${params.trackpy_percentile} -trackpy_search_range 9 -peak_separation ${params.trackpy_separation} -stem ${stem}
    """
}

process Process_tracks {
    echo true
    storeDir params.out_dir + "/anchor_peaks"
    /*publishDir params.out_dir + "/anchor_peaks", mode:"copy"*/

    input:
    tuple val(stem), file(tsv)

    output:
    tuple val(stem), file("${stem}_processed_peaks*"), emit: processed_tracks
    /*file("*png"), optional*/

    script:
    """
    python ${workflow.projectDir}/py_scripts/process_tracks.py -tracks ${tsv} -stem ${stem}
    """
}

process Extract_peak_intensities_from_preprocessed_arrays {
    echo true
    storeDir params.out_dir + "peak_intensities"
    /*publishDir params.out_dir + "peak_intensities", mode:"copy"*/

    input:
    tuple val(stem), file(peaks)
    file(imgs)

    output:
    tuple val(stem), file("${stem}_extracted_peak_intensities*"), file("spot_locs.csv"), emit: paeks_for_decoding

    script:
    """
    ls
    #python ${workflow.projectDir}/py_scripts/extract_peak_intensities.py -zarr ${imgs} -anchors ${peaks} -out ${stem}
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
    /*println Get_meatdata.out[0]*/
    extract_anchor_chs(channel.fromPath(params.ome_tif), params.anchor_ch_indexes, params.format)
    Enhance_spots(extract_anchor_chs.out.extracted_chs)
    Call_peaks_in_anchor(Enhance_spots.out.processed_anchor_img_zarr)
    /*println Call_peaks_in_anchor.out.peaks_from_anchor_ch*/
    Process_tracks(Call_peaks_in_anchor.out.peaks_from_anchor_ch)
    Extract_peak_intensities_from_preprocessed_arrays(Process_tracks.out,
        channel.fromPath(params.ome_tif))
}
