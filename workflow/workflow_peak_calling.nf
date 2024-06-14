#!/usr/bin/env/ nextflow

process Spotiflow_call_peaks {
    cache true

    label "gpu_normal"

    container 'bioinfotongli/spotiflow:latest'
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    storeDir params.out_dir

    input:
    tuple val(meta), path(root)

    output:
    tuple val(meta), path("peaks_Y*_X*.csv"), emit: peaks

    script:
    def args = task.ext.args ?: ''
    """
    Spotiflow_call_peaks.py run \
        -image_path ${root}/0 \
        ${args}
    """
}


process Spotiflow_merge_peaks {
    debug true
    cache true

    container 'bioinfotongli/spotiflow:latest'
    storeDir params.out_dir

    input:
    tuple val(meta), path(csvs)

    output:
    tuple val(meta), path("merged_peaks.wkt"), emit: merged_peaks

    script:
    def args = task.ext.args ?: ''
    """
    Spotiflow_post_process.py run \
        ${args} \
        ${csvs} \
    """
}

workflow Spotiflow_run {
    take:
    zarrs

    main:
    Spotiflow_call_peaks(zarrs)
    Spotiflow_merge_peaks(Spotiflow_call_peaks.out.peaks.collect())
}