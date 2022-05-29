#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process PREFIX_FILENAME {

    // maxForks 1

    // echo true

    input:
        path(fasta)
        val(prefix_fn)

    output:
        tuple stdout,
            path(fasta)

    script:
    """
    prefix_filename.py ${fasta} ${prefix_fn} 
    """
}


process REFORMAT_FASTA {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            path(fasta)

    output:
        tuple val(prefix_filename),
            stdout,
            path("reformatted.fa")

    script:
    """
    reformat_fasta.py ${fasta}
    """
}


process FASTA_VARIANT_SITES {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("reformatted.fa")

    output:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("reformatted.fa"),
            path("ref_variants_in_fasta.csv"),
            path("variants_in_fasta.csv")
        
    script:
    """
    fasta_variant_sites.py reformatted.fa
    """
}


process PAIRWISE_TABLE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("reformatted.fa"),
            path("ref_variants_in_fasta.csv"),
            path("variants_in_fasta.csv")

    output:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("variants_in_fasta.csv"),
            path("pairwise_table.csv")

    script:
    """
    pairwise_table_fasta.py variants_in_fasta.csv reformatted.fa
    """
}


process PAIRWISE_BIALLELIC_TABLE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("variants_in_fasta.csv"),
            path("pairwise_table.csv")

    output:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("variants_in_fasta.csv"),
            path("pairwise_biallelic_table.csv")

    script:
    // NOTE
    // There was serious error which prevented it from working properly, but is now fixed.
    // Need to fix in other projects that use this also
    """
    biallelic_filter.py pairwise_table.csv
    """
}

process WATTERSON_ESTIMATE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size_seq_len),
            path("variants_in_fasta.csv"),
            path("pairwise_biallelic_table.csv")

    output:
        tuple val(prefix_filename),
            stdout,
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv")

    script:
    """
    snps=\$(wc -l variants_in_fasta.csv | awk '{split(\$0,a," "); print a[1]}')
    sample_size=\$(echo ${sample_size_seq_len} | cut -d',' -f1)
    genome_len=\$(echo ${sample_size_seq_len} | cut -d',' -f2)
    original_watterson.py \$genome_len \$snps \$sample_size
    """
}


process LOOKUP_TABLE_LDPOP {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv"),
            path("lookupTable.txt")

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    sample_size=\$(echo ${sample_size_seq_len} | cut -d',' -f1)
    ldtable.py --cores $task.cpus -n \$sample_size -th ${theta_est} -rh ${params.lookup_grid} --approx > lookupTable.txt
    """
}


process DOWNSAMPLED_LOOKUP_TABLE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv")
        
        val theta
        path downsampled_lookup_tables
        

    output:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv"),
            path("lookupTable.txt")

    script:
    """
    sample_size=\$(echo ${sample_size_seq_len} | cut -d',' -f1)
    reformat_downsampled_lk_table.py lk_downsampled_\$sample_size.csv \$sample_size ${theta} ${params.lookup_grid}
    """
}


process PAIRWISE_LOOKUP_FORMAT {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv"),
            path("lookupTable.txt")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv"),
            path("lookupTable.txt"),
            path("lookup_format.csv")

    script:
    """
    pairwise_lookup_format_pyrho.py pairwise_biallelic_table.csv
    """
}


process CUSTOM_HAP_SETS_AND_MERGE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("pairwise_biallelic_table.csv"),
            path("lookupTable.txt"),
            path("lookup_format.csv")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("lookupTable.txt"),
            path("table_ids_for_eq3.csv"),
            path("eq3.csv")

    script:
    """
    custom_hap_sets_and_merge.py lookupTable.txt pairwise_biallelic_table.csv lookup_format.csv ${params.lookup_grid} > table_ids_for_eq3.csv
    """
}


process P_IJ_GRID {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size_seq_len),
            path("lookupTable.txt"),
            path("table_ids_for_eq3.csv"),
            path("eq3.csv")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            path("lookupTable.txt"),
            path("table_ids_for_eq3.csv"),
            path("eq3.csv"),
            path("p_ij_grid.csv")

    script:
    """
    genome_len=\$(echo ${sample_size_seq_len} | cut -d',' -f2)
    pij_grid_vectorised.py \$genome_len ${params.recom_tract_len} ${params.lookup_grid} eq3.csv
    """
}


process PAIRWISE_ESTIMATOR {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}
    
    input:
        tuple val(prefix_filename),
            val(theta_est),
            path("lookupTable.txt"),
            path("table_ids_for_eq3.csv"),
            path("eq3.csv"),
            path("p_ij_grid.csv")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            path("collected_likelihoods.csv")
    
    script:
    """
    pairwise_rho_estimator_intp_rect_biv.py eq3.csv table_ids_for_eq3.csv p_ij_grid.csv lookupTable.txt ${params.lookup_grid}
    """

}


process FINAL_RESULTS {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            path("collected_likelihoods.csv")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            path("final_results.txt")

    script:
    """
    final_results.py collected_likelihoods.csv ${theta_est}
    """
}


process PROCESS_OUTPUT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            path("final_results.txt")

    output:
        path "processed_results.csv", emit: processed_results_csv
        path "theta_est.csv", emit: theta_est_csv

    script:
        """
        custom_est_process_output.py final_results.txt
        echo ${theta_est} > theta_est.csv
        """

}

process PLOT_RESULTS{
    publishDir "Output/Results", mode: "copy"

    maxForks 1

    input:
        path collectedFile


    output:
        path "rho_comparision.png", emit: rho_comparision_png
        path "max_lk_comparision.png", emit: max_lk_comparision_png

    script:
        """
        plot_results.py collected_results.csv
        """

}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // params.theta = 0.01 // scaled theta
    params.recom_tract_len = 1000
    params.lookup_grid = "101,100" // The range of rho values used to generate lookup tables

    params.prefix_filename = 'none' // prefix string to output filenames to help distinguish runs
    params.input_fasta = 'none'
    // params.lookup_tables = "Lookup_tables"
    // params.lookup_tables = "/Volumes/Backup/Lookup_tables/Lookup_tables_m_0.01_r_0-100"
    // params.lookup_tables = "/shared/homes/11849395/Lookup_tables/Lookup_tables_0-100"

    // Input verification
    if (params.input_fasta == 'none') {
        println "No input .fa specified. Use --input_fasta [.fa]"
        exit 1
    }

    // Channels
    input_fasta_channel = Channel.fromPath( params.input_fasta )
    // downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv" ).collect()

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need
    PREFIX_FILENAME(input_fasta_channel, params.prefix_filename)

    REFORMAT_FASTA(PREFIX_FILENAME.out)

    FASTA_VARIANT_SITES(REFORMAT_FASTA.out)

    PAIRWISE_TABLE(FASTA_VARIANT_SITES.out)

    PAIRWISE_BIALLELIC_TABLE(PAIRWISE_TABLE.out)

    WATTERSON_ESTIMATE(PAIRWISE_BIALLELIC_TABLE.out)

    LOOKUP_TABLE_LDPOP(WATTERSON_ESTIMATE.out)

    // DOWNSAMPLED_LOOKUP_TABLE(WATTERSON_ESTIMATE.out, params.theta, downsampled_lookup_tables)

    PAIRWISE_LOOKUP_FORMAT(LOOKUP_TABLE_LDPOP.out)

    CUSTOM_HAP_SETS_AND_MERGE(PAIRWISE_LOOKUP_FORMAT.out)

    P_IJ_GRID(CUSTOM_HAP_SETS_AND_MERGE.out)

    PAIRWISE_ESTIMATOR(P_IJ_GRID.out)

    FINAL_RESULTS(PAIRWISE_ESTIMATOR.out)

    PROCESS_OUTPUT(FINAL_RESULTS.out)

    // collectedFile = PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    // PLOT_RESULTS(collectedFile)
}