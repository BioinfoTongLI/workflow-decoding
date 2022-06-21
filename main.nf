#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

// minimal parameter set
params.ome_tif = ''
params.out_dir = ''
params.rna_spot_size = [5, 7]
params.anchor_ch_indexes = 2

// needed when decoding
params.taglist_name = "taglist.csv"
params.channel_info_name = "channel_info.csv"
params.codebook = ""

params.anchor_peaks_tsv = "" // if avaiable peaks were detected by Synquant

params.max_n_worker = 28

params.trackpy_percentile = [85, 90]
params.trackpy_separation = 2
params.trackpy_search_range = 5


// not used in this version
/*params.tile_name : "N1234F_tile_names.csv"*/
/*params.coding_ch_starts_from = 0*/
/*params.anchor_available = 1*/

/*params.known_anchor = "c01 Alexa 647"*/

/*
 * bf2raw: The bioformats2raw application converts the input image file to
 * an intermediate directory of tiles in the output directory.
 */
process bf2raw {
    debug true
    container "openmicroscopy/bioformats2raw:0.4.0"
    storeDir params.out_dir + "/raws"
    /*publishDir params.out_dir, mode:"copy"*/

    input:
    path(img)

    output:
    tuple val(stem), file("${stem}")

    script:
    stem = img.baseName
    """
    /opt/bioformats2raw/bin/bioformats2raw --max_workers ${params.max_n_worker} --no-hcs $img "${stem}"
    """
}


process Codebook_conversion {
    echo true
    cache "lenient"
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    storeDir params.out_dir + "/decoding_metadata"

    input:
    file(codebook)

    output:
    path "taglist.csv", emit: taglist_name
    path "channel_info.csv", emit: channel_info_name
    /*path "channel_info.pickle", emit: channel_infos*/

    script:
    """
    codebook_convert.py -csv_file ${codebook}
    """
}


process Get_meatdata {
    echo true
    cache "lenient"
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    storeDir params.out_dir + "/decoding_metadata"
    /*publishDir params.out_dir + "/decoding_metadata", mode:"copy"*/

    input:
    path(taglist_name)
    path(channel_info_name)

    output:
    path "barcodes_01.npy", emit: barcodes
    path "gene_names.npy", emit: gene_names
    path "channel_info.pickle", emit: channel_infos

    script:
    """
    get_metadata.py -auxillary_file_dir ./  -taglist_name ${taglist_name} -channel_info_name ${channel_info_name}
    """
}


process Enhance_spots {
    debug true
    cache "lenient"
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    containerOptions "--gpus all"
    storeDir params.out_dir + "/enhanced_anchor_channels"
    /*publishDir params.out_dir + "/anchor_spots", mode:"copy"*/

    input:
    tuple val(stem), path(zarr)
    val anchor_ch_indexes
    each rna_spot_size

    output:
    tuple val(rna_spot_size), val(stem), path("${stem}_spot_enhanced_diam_${rna_spot_size}"), emit:ch_with_peak_img

    script:
    """
    helper.py enhance_all --diam ${rna_spot_size} --zarr_in ${zarr}/0 --stem ${stem}
    """
}


process Deepblink_and_Track {
    echo true
    cache "lenient"
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:deepblink"
    containerOptions "--gpus all -v ${workflow.projectDir}:${workflow.projectDir}"
    /*storeDir params.out_dir + "/anchor_spots"*/
    publishDir params.out_dir + "/anchor_spots", mode:"copy"

    input:
    tuple val(stem), path(zarr)

    output:
    tuple val(stem), file("${stem}_max_*_peaks.tsv"), emit: peaks_from_anchor_chs

    script:
    """
    python3 ${workflow.projectDir}/py_scripts/deepblink_wrap.py --zarr_in ${zarr}/0 --stem ${stem} --tpy_search_range 3
    """
}


process Call_peaks_in_anchor {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    containerOptions "--gpus all"
    storeDir params.out_dir + "/anchor_spots"
    /*publishDir params.out_dir + "/anchor_spots", mode:"copy"*/

    maxForks 1

    input:
    tuple val(rna_spot_size), val(stem), file(anchor_zarr)
    val(anchor_ch_index)
    each percentile
    each separation
    each search_range

    output:
    tuple val(rna_spot_size), path("${stem}_detected_peaks_diam_${rna_spot_size}_percentile_${percentile}_sep_${separation}_search_range_${search_range}.tsv"), emit: peaks_from_anchor_chs
    /*tuple val(stem), file("${stem}_tracked_peaks.tsv"), emit: tracked_peaks_from_anchor_chs*/

    script:
    """
    helper.py call_peaks --zarr_in ${anchor_zarr}/0/${anchor_ch_index} --stem ${stem} --diam ${rna_spot_size} --tp_percentile ${percentile} --peak_separation ${separation} --tpy_search_range ${search_range} --anchor_ch_index ${anchor_ch_index}
    """
}


process Process_peaks {
    echo true
    cache "lenient"
    container "docker://rapidsai/rapidsai:cuda11.2-base-ubuntu20.04-py3.8"
    containerOptions "--gpus all -v ${workflow.projectDir}:${workflow.projectDir}"
    /*storeDir params.out_dir + "/anchor_spots"*/
    publishDir params.out_dir + "/anchor_spots", mode:"copy"

    input:
    tuple val(stem), path(track_tsv)

    output:
    tuple val(stem), file("${stem}_processed_tracks.tsv"), emit: tracked_peaks

    script:
    """
    /opt/conda/envs/rapids/bin/python ${workflow.projectDir}/py_scripts/process_tracks.py --tsv ${track_tsv} --stem ${stem}
    """
}


process Extract_peak_intensities {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    containerOptions "--gpus all"
    /*storeDir params.out_dir + "/peak_intensities"*/
    publishDir params.out_dir + "/peak_intensities", mode:"copy"

    input:
    tuple val(rna_spot_size), path(peaks), val(stem), path(imgs), path(channel_info)

    output:
    tuple val(new_stem), file("${new_stem}_extracted_peak_intensities.npy"), file("${new_stem}_peak_locs.csv"), emit: peaks_for_preprocessing

    script:
    new_stem = peaks.baseName
    """
    extract_peak_intensities.py --raw_zarr ${imgs}/0 \
        --peaks ${peaks} \
        --stem ${new_stem} \
        --channel_info ${channel_info} \
        --coding_cyc_starts_from 1
    """
}


process Preprocess_peak_profiles {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    containerOptions "--gpus all"
    storeDir params.out_dir + "/preprocessed_peak_intensities"
    /*publishDir params.out_dir + "/preprocessed_peak_intensities", mode:"copy"*/

    input:
    tuple val(stem), file(profiles), file(peak_locations)

    output:
    tuple val(stem), file("${stem}_filtaered_peak_intensities.npy"), file("${stem}_filtered_peak_locs.csv"), emit: peaks_for_decoding

    script:
    """
    preprocess_peak_profiles.py --profiles ${profiles} --stem ${stem} --spot_loc ${peak_locations}
    """
}


process Decode_peaks {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    containerOptions "--gpus all"
    storeDir params.out_dir + "decoded"
    /*publishDir params.out_dir + "/decoded", mode:"copy"*/

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
    decode.py --spot_profile ${spot_profile} --spot_loc ${spot_loc} --barcodes_01 ${barcodes_f} --gene_names ${gene_names_f} --channels_info ${channel_info_f} --stem ${stem}
    """
}


process Heatmap_plot {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest"
    publishDir params.out_dir + "/plots", mode:"copy"

    input:
    path(decoded_df)
    /*file decode_out_parameters_f from decode_out_parameters*/
    /*file channel_info_f from channel_info_for_plot*/

    output:
    path("*.png")

    script:
    """
    do_plots.py --decoded_df $decoded_df
    #-decode_out_params decode_out_parameters_f -channels_info {channel_info_f}
    """
}


workflow {
    Codebook_conversion(Channel.fromPath(params.codebook))
    Get_meatdata(Codebook_conversion.out.taglist_name, Codebook_conversion.out.channel_info_name)
    peak_calling()
    Extract_peak_intensities(
        peak_calling.out.peaks_and_enhanced_raw.combine(Get_meatdata.out.channel_infos)
    )
    Preprocess_peak_profiles(Extract_peak_intensities.out.peaks_for_preprocessing)
    Decode_peaks(
        Preprocess_peak_profiles.out.peaks_for_decoding,
        Get_meatdata.out.barcodes,
        Get_meatdata.out.gene_names,
        Get_meatdata.out.channel_infos
    )
    /*Heatmap_plot(Decode_peaks.out[0])*/
}

workflow peak_calling {
    bf2raw(params.ome_tif)
    Enhance_spots(bf2raw.out, params.anchor_ch_indexes, channel.from(params.rna_spot_size))
    if (params.anchor_peaks_tsv != "") {
        peaks = Channel.fromPath(params.anchor_peaks_tsv)
    } else {
        Call_peaks_in_anchor(Enhance_spots.out.ch_with_peak_img, params.anchor_ch_indexes,
            channel.from(params.trackpy_percentile),
            channel.from(params.trackpy_separation),
            channel.from(params.trackpy_search_range)
        )
        peaks = Call_peaks_in_anchor.out.peaks_from_anchor_chs
        /*Process_peaks(Deepblink_and_Track.out.peaks_from_anchor_chs)*/
    }
    Enhance_spots.out.ch_with_peak_img
        .flatMap{it-> [it] * params.trackpy_percentile.size()}
        .set{duplicated_spots}
        /*.view()*/
    peaks.cross(duplicated_spots)
        .map{it -> [it[0][0], it[0][1], it[1][1], it[1][2]]}
        .set{peaks_and_enhanced_raw}
        /*.view()*/

    emit:
    peaks_and_enhanced_raw = peaks_and_enhanced_raw

}
