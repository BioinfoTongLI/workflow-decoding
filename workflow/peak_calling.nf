process Deepblink_and_Track {
    debug true
    cache "lenient"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'singularity-image':
        'docker-image'}"
    containerOptions "${ workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"

    storeDir params.out_dir + "/anchor_spots"

    input:
    tuple val(stem), path(zarr)

    output:
    tuple val(stem), file("${stem}_max_*_peaks.tsv"), emit: peaks_from_anchor_chs

    script:
    """
    deepblink_wrap.py --zarr_in ${zarr}/0 --stem ${stem} --tp_search_range 3
    """
}

process Process_tracks {
    debug true
    cache "lenient"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'singularity-image':
        'docker-image'}"
    containerOptions "${ workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"

    storeDir params.out_dir + "/anchor_spots"

    input:
    tuple val(stem), path(track_tsv)

    output:
    tuple val(stem), file("${stem}_processed_tracks.tsv"), emit: tracked_peaks

    script:
    """
    process_tracks.py --tsv ${track_tsv} --stem ${stem}
    """
}
