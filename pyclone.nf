#!/usr/bin/env nextflow

process pyclone_fit {
    tag "$case_id"
    conda '/home/data/data_dingyangliu/miniconda3/envs/pyclone-vi'

    publishDir "pyclone_output/${case_id}/", mode: 'copy'

    input:
    tuple val(case_id), path(input_file)

    output:
    tuple val(case_id), path("${case_id}.h5")

    script:
    """
    pyclone-vi fit \
        -i ${input_file} \
        -o ${case_id}.h5 \
        -c 40 \
        -d beta-binomial \
        -r 10
    """
}

process pyclone_write {
    tag "$case_id"
    conda '/home/data/data_dingyangliu/miniconda3/envs/pyclone-vi'
    publishDir "pyclone_output/${case_id}/", mode: 'copy'

    input:
    tuple val(case_id), path(input_h5)

    output:
    path("${case_id}.pyclone.tsv")

    script:
    """
    pyclone-vi write-results-file \
        -i ${input_h5} \
        -o ${case_id}.pyclone.tsv
    """
}

workflow {
    Channel
    //修改fromPath路徑
        .fromPath("data/*.tsv")
        .map { file -> 
            def case_id = file.getBaseName().replaceFirst(/\.tsv$/, "")
            tuple(case_id, file)
        }
        .set { input_samples }
    view(input_samples)
    fit_result = pyclone_fit(input_samples)
    pyclone_write(fit_result)
}
