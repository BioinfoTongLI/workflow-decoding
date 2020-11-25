#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

params.auxillary_file_dir = "/nfs/team283_imaging/NT_ISS/playground_Tong/KR0018/new_opt/gmm-input/"
params.tile_name = "tile_names.csv"
params.taglist_name = "taglist.csv"
params.channel_info_name = "channel_info.csv"
params.data_dir = '/nfs/team283_imaging/NT_ISS/playground_Tong/KR0018/new_opt/tiles/splits/selected_tiles/'
params.proj_ID = "mouse_brain"
params.out_dir = "./"
params.tile_size = 1000
params.filename_prefix = "out_opt_flow_registered"
params.tile_name_pattern = "{prefix}_{tile_xy}_c0{cycle_i}_{ch_name}.tif"
params.anchor_available = 1
params.trackpy_percentile = 64
params.coding_ch_starts_from = 0


process get_meatdata {
    echo true
    cache "lenient"
    containerOptions = "-B /nfs:/nfs:ro"

    output:
    path "barcodes_01.npy" into barcodes_01
    path "gene_names.npy" into gene_names
    path "channel_info.pickle" into channels_info, channel_info_for_decode

    script:
    """
    python ${workflow.projectDir}/get_metadata.py -auxillary_file_dir $params.auxillary_file_dir  -taglist_name ${params.taglist_name} -channel_info_name ${params.channel_info_name}
    """
}


process get_spots {
    echo true
    cache "lenient"
    storeDir params.out_dir + "/" + params.proj_ID + "_trackpy_spots"
    /*publishDir params.out_dir + "/trackpy_spots", mode:"copy"*/

    input:
    path channel_info_f from channels_info

    output:
    path "*_spots_trackpy_spot_profile.npy" into spot_profile
    path "*_spots_trackpy_locations.csv" into spot_locations

    script:
    """
    python ${workflow.projectDir}/get_spots.py -tifs_dir $params.data_dir -stem ${params.proj_ID} -tile_name ${params.auxillary_file_dir}/${params.tile_name} -tile_size ${params.tile_size} -filename_prefix ${params.filename_prefix} -tile_name_pattern ${params.tile_name_pattern} -anchor_available ${params.anchor_available} -trackpy_percentile ${params.trackpy_percentile} -channel_info ${channel_info_f} -coding_ch_starts_from ${params.coding_ch_starts_from}
    """
}


process decode {
    echo true
    /*storeDir params.out_dir + "/" + params.proj_ID + "_decoded"*/
    publishDir params.out_dir + "/" + params.proj_ID + "_decoded", mode:"copy"

    input:
    path spot_profile from spot_profile
    path spot_loc from spot_locations
    path barcodes_f from barcodes_01
    path gene_names_f from gene_names
    path channel_info_f from channel_info_for_decode

    output:
    path "${params.proj_ID}_decoded_df.tsv" into decoded_df
    path "${params.proj_ID}_decode_out_parameters.pickle" into decode_out_parameters

    script:
    """
    python ${workflow.projectDir}/decode.py -spot_profile ${spot_profile} -spot_loc ${spot_loc} -barcodes_01 ${barcodes_f} -gene_names ${gene_names_f} -channels_info ${channel_info_f} -stem ${params.proj_ID}
    """
}
