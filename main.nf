params.ref = "/home/data/data_thousand/gatk_index/Homo_sapiens_assembly38.fasta"
params.gc  = "/home/data/data_dingyangliu/ASC_0217/test_sequenza/genome_gc50.wig.gz"


process bam2seqz {
    tag "$sample_id"
    conda '/home/data/data_dingyangliu/miniconda3/envs/seqz'
    publishDir "output/${sample_id}", mode: 'copy' 

    input:  
    tuple val(sample_id), path(tumor_bam), path(normal_bam)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.seqz.gz")

    script:
    """
    mkdir -p ${sample_id}
    sequenza-utils bam2seqz \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -gc ${params.gc} \\
        -F ${params.ref} \\
        -o ${sample_id}/${sample_id}.seqz.gz
    """
}

process bin {
    tag "$sample_id"
    conda '/home/data/data_dingyangliu/miniconda3/envs/seqz'
    publishDir "output/${sample_id}", mode: 'copy' 

    input:  
    tuple val(sample_id), path(seqz)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_bin50.seqz.gz")

    script:
    """
    mkdir -p ${sample_id}
    sequenza-utils seqz_binning --seqz $seqz -w 50 -o ${sample_id}/${sample_id}_bin50.seqz.gz
    """
    
}

process scarHRD {
    tag "$sample_id"
    conda '/home/data/data_dingyangliu/miniconda3/envs/scarHRD'
    publishDir "output/${sample_id}", mode: 'copy' 

    input:
    tuple val(sample_id), path(seqz_file)

    output:
    tuple val(sample_id), path("scarHRD_result")

    script:
    """
    Rscript ${projectDir}/bin/run_scarHRD.R ${sample_id} ${seqz_file}
    """
}

process sequenza {
    tag "$sample_id"
    conda '/home/data/data_dingyangliu/miniconda3/envs/r-sequenza'
    publishDir "output/${sample_id}", mode: 'copy' 
    cpus 20

    input:
    tuple val(sample_id), path(seqz_file)

    output:
    tuple val(sample_id), path("sequenza_result")

    script:
    """
    Rscript ${projectDir}/bin/run_sequenza.R ${sample_id} ${seqz_file} ${task.cpus}
    """
    }

process merge_hrd_results {
    tag "merge HRD results"
    publishDir "output", mode: 'copy' 

    input:
    path hrd_files

    output:
    path "HRD_summary.tsv"

    script:
    """
    echo -e "Sample\\tHRD\\tTelomeric_AI\\tLST\\tHRD_sum" > HRD_summary.tsv
    for f in ${hrd_files}; do
        sample_name=\$(basename "\$f" "_bin50.seqz._HRDresults.txt")
        tail -n 1 "\$f" | awk -v s="\$sample_name" '{print s "\\t" \$0}' >> HRD_summary.tsv
    done
    """
    }

process merge_sequenza_results {
    tag "merge Sequenza results"
    publishDir "output", mode: 'copy' 

    input:
    path sequenza_files

    output:
    path("all_segments_merged.tsv")

    script:
    """
    echo -e "chromosome\tstart.pos\tend.pos\tBf\tN.BAF\tsd.BAF\tdepth.ratio\tN.ratio\tsd.ratio\tCNt\tA\tB\tLPP\tcellularity\tploidy\tsample.id" > all_segments_merged.tsv
    for f in ${sequenza_files}; do
        tail -n +2 \$f >> all_segments_merged.tsv
    done
    """
}



// ✅ Pipeline 結束時發送通知
workflow.onComplete {
    def status = workflow.success ? "✅ 成功" : "❌ 失敗"
    def summary = "【Sequenza 工作流程完成】\n狀態: ${status}\n耗時: ${workflow.duration}"
    println summary
    ["python3", "send_msg.py", summary.replaceAll('\n', ' ')].execute().waitFor()
}

workflow {
    def msg = "【Sequenza 工作流程啟動】\n時間: ${new Date()}"
    println msg
    ["python3", "send_msg.py", msg.replaceAll('\n', ' ')].execute().waitFor()


    Channel
        .fromPath("tn_pairs.tsv")
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.tumor_bam), file(row.normal_bam)) }
        //.take(1) // 如果只測試一個
        .set { tumor_normal_pairs }

    bam2seqz(tumor_normal_pairs)
    bin(bam2seqz.out)
    sequenza(bin.out)
    scarHRD(bin.out)

    // merge HRD results
    scarHRD.out
    .map { sample_id, result_dir -> 
        file("${result_dir}/${sample_id}_bin50.seqz._HRDresults.txt")
    }
    .collect()
    .set { all_hrd_results }

    merge_hrd_results(all_hrd_results)

    // merge Sequenza results
    sequenza.out
    .map { sample_id, result_dir ->         
        file("${result_dir}/${sample_id}_segments_modified.txt")
    }
    .collect()
    .set { all_sequenza_results }

    merge_sequenza_results(all_sequenza_results)
}
