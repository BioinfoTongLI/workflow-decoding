#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tif = 'path/to/ome.tiff'
params.out_dir = "./test/"
params.known_anchor = "c01 Alexa 647"
params.trackpy_separation = 2
params.rna_spot_size = 5
params.trackpy_percentile = 90
params.anchor_ch_indexes = 4
params.format = "--zarr"

params.decode = true
params.auxillary_file_dir = "/nfs/team283_imaging/NT_ISS/playground_Tong/KR0018/new_opt/gmm-input/"
params.taglist_name = "taglist.csv"
params.channel_info_name = "channel_info.csv"

// not used in this version
params.coding_ch_starts_from = 0
params.anchor_available = 1
params.max_n_worker = 25

/*ome_tif_ch = Channel.fromPath(params.ome_tif).*/
    /*into{ome_tif_for_anchor_peak_calling; ome_tif_for_peak_intensity_extraction}*/

/*
 * bf2raw: The bioformats2raw application converts the input image file to
 * an intermediate directory of tiles in the output directory.
 */
process bf2raw {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/bayraktar-lab/image-convert:latest"
    /*storeDir params.out_dir + "/raws"*/
    publishDir params.out_dir, mode:"copy"

    input:
    path(img)

    output:
    tuple val(stem), file("${stem}")

    script:
    stem = img.baseName
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --no-hcs $img "${stem}"
    """
}


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


process Enhance_spots {
    echo true
    cache "lenient"
    storeDir params.out_dir + "/anchor_spots"
    /*publishDir params.out_dir + "/anchor_spots", mode:"copy"*/
    containerOptions "--nv"

    input:
    tuple val(stem), path(zarr)
    val anchor_ch_indexes
    path channel_info

    output:
    tuple val(stem), path("${stem}_spot_enhanced"), emit:ch_with_peak_img
    /*tuple val(stem), path("${stem}_spot_enhanced.tif"), emit:ch_with_peak_tif*/

    script:
    """
    python ${workflow.projectDir}/helper.py enhance_spots --diam ${params.rna_spot_size} --ch_info ${channel_info}  --zarr_in ${zarr}/0 --stem ${stem} --anchor_ch_ind ${anchor_ch_indexes}
    """
}


process Deepblink_and_Track {
    echo true
    cache "lenient"
    container "/home/ubuntu/sifs/deepblink-2021-07-22-caf58c80de23.sif"
    containerOptions "--nv"
    storeDir params.out_dir + "/anchor_spots"
    /*publishDir params.out_dir + "/anchor_spots", mode:"copy"*/

    input:
    tuple val(stem), path(zarr)

    output:
    tuple val(stem), file("${stem}_max_*_peaks.tsv"), emit: peaks_from_anchor_chs

    script:
    """
    python3 ${workflow.projectDir}/deepblink_wrap.py --zarr_in ${zarr}/0 --stem ${stem} --tpy_search_range 3
    """
}


process Call_peaks_in_anchor {
    echo true
    /*storeDir params.out_dir + "/anchor_spots"*/
    publishDir params.out_dir + "/anchor_spots", mode:"copy"

    input:
    tuple val(stem), file(anchor_zarr)

    output:
    tuple val(stem), file("${stem}_tracked_peaks.tsv"), emit: peaks_from_anchor_chs

    script:
    """
    python ${workflow.projectDir}/helper.py call_peaks --zarr_in ${anchor_zarr}/0 --stem ${stem} --diam ${params.rna_spot_size} --tp_percentile ${params.trackpy_percentile} --peak_separation ${params.trackpy_separation} --tpy_search_range 5
    # serach range fixed to 5 as tp could not handle more depth
    """
}


process Process_peaks {
    echo true
    cache "lenient"
    /*container "container/rapids_gmm.sif"*/
    container "docker://rapidsai/rapidsai:cuda11.2-base-ubuntu20.04-py3.8"
    containerOptions "--nv"
    /*storeDir params.out_dir + "/anchor_spots"*/
    publishDir params.out_dir + "/anchor_spots", mode:"copy"

    input:
    tuple val(stem), path(track_tsv)

    output:
    tuple val(stem), file("${stem}_processed_tracks.tsv"), emit: tracked_peaks

    script:
    """
    /opt/conda/envs/rapids/bin/python ${workflow.projectDir}/process_tracks.py --tsv ${track_tsv} --stem ${stem}
    """
}


process Extract_peak_intensities {
    echo true
    /*storeDir params.out_dir + "/peak_intensities"*/
    publishDir params.out_dir + "/peak_intensities", mode:"copy"

    input:
    tuple val(stem), file(peaks), file(imgs)
    file channel_info

    output:
    tuple val(stem), file("${stem}_extracted_peak_intensities*"), file("${stem}_peak_locs.csv"), emit: peaks_for_decoding

    script:
    """
    python ${workflow.projectDir}/extract_peak_intensities.py --raw_zarr ${imgs}/0 --peaks ${peaks} --stem ${stem} --channel_info ${channel_info} --coding_cyc_starts_from 1
    """
}


process Decode_peaks {
    echo true
    /*storeDir params.out_dir + "decoded"*/
    publishDir params.out_dir + "decoded", mode:"copy"
    containerOptions "--nv"

    when:
    params.decode

    input:
    tuple val(stem), file(spot_profile), file(spot_loc)
    file barcodes_f
    file gene_names_f
    file channel_info_f

    output:
    file "${stem}_decoded_df.tsv"
    file "${stem}_decode_out_parameters.pickle"

    script:
    """
    python ${workflow.projectDir}/decode.py --spot_profile ${spot_profile} --spot_loc ${spot_loc} --barcodes_01 ${barcodes_f} --gene_names ${gene_names_f} --channels_info ${channel_info_f} --stem ${stem}
    """
}


process Do_Plots {
    echo true
    containerOptions " -v " + baseDir + ":/gmm_decoding/:ro"
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
    bf2raw(params.ome_tif)
    Enhance_spots(bf2raw.out, params.anchor_ch_indexes, Get_meatdata.out.channel_infos)
    Deepblink_and_Track(Enhance_spots.out.ch_with_peak_img)
    /*Call_peaks_in_anchor(Enhance_spots.out.ch_with_peak_img)*/
    /*Process_peaks(Deepblink_and_Track.out.peaks_from_anchor_chs)*/
    Extract_peak_intensities(
        Deepblink_and_Track.out.peaks_from_anchor_chs.join(bf2raw.out),
        Get_meatdata.out.channel_infos
    )
    Decode_peaks(
        Extract_peak_intensities.out.peaks_for_decoding,
        Get_meatdata.out.barcodes,
        Get_meatdata.out.gene_names,
        Get_meatdata.out.channel_infos
    )
}
