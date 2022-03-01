#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    maxForks 1 // Run sequentially

    input:
        val rho_rate
        val theta
        val genome_size
        val sample_size
        val seed

    output:
        val "${rho_rate}", emit: rho_rate
        val "${theta}", emit: theta
        val "${genome_size}", emit: genome_size
        val "${sample_size}", emit: sample_size
        val "${seed}", emit: seed
        
        val "rho_${rho_rate}_theta_${theta}_genome_size_${genome_size}_sample_size_${sample_size}_seed_${seed}", emit: path_fn_modifier

    script:
    """
    """

}


process MS {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 

    input:
        val rho_rate
        val sample_size
        val seed
        val genome_size
        val path_fn_modifier

    output:
        path "trees.txt", emit: trees_txt
  
    script:
    """
    ms ${sample_size} 1 -T -seeds ${seed}  -t ${params.mutation_rate} -r ${rho_rate} ${genome_size} -c 10 ${params.recom_tract_len} > trees.txt
    """
}


process FAST_SIM_BAC {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val theta
        val seed
        val path_fn_modifier

    output:
        // path "rho_calc.txt", emit: rho_rho_calc_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${seed} -T -t ${theta} -r ${rho_rate} ${params.recom_tract_len} > trees.txt
    """
}


process MS_PRIME {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "rho_calculation.txt", emit: rho_rho_calculation_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    run_msprime.py
    calc_rho.py ${params.effective_pop_size} ${rho_rate} ${genome_size}
    """
}


process CLEAN_TREES {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}
    
    maxForks 1

    input:
        path trees
        val path_fn_modifier


    output:
        path "cleanTrees.txt", emit: cleanTrees_txt

    script:
    """
    clean_trees.py trees.txt
    """
}


process SEQ_GEN {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path cleanTrees
        val genome_size
        val seed
        val path_fn_modifier

    output:
        path "seqgenOut.fa", emit: seqgenout_fa

    script:
    // 1 partiion per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


// process ALLELE_FREQ {
//     publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

//     maxForks 1

//     input:
//         path seqgenOut
//         val path_fn_modifier

//     output:
//         path "allele_frequencies.txt", emit: allele_frequencies_txt

//     script:
//     """
//     allele_frequencies.py seqgenOut.fa
//     """
// }


process REFORMAT_FASTA {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path seqgenOut
        val path_fn_modifier

    output:
        path "reformatted.fa", emit: reformatted_fa

    script:
    """
    reformat_fasta.py seqgenOut.fa
    """
}


process FASTA_VARIANT_SITES {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path reformatted_fa
        val path_fn_modifier

    output:
        path "ref_variants_in_fasta.csv", emit: ref_variants_in_fasta_csv
        path "variants_in_fasta.csv", emit: variants_in_fasta_csv
        
    script:
    """
    fasta_variant_sites.py reformatted.fa
    """
}


process PAIRWISE_TABLE {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path variants_in_fasta_csv
        path reformatted_fa
        val single_end_restrict
        val path_fn_modifier

    output:
        path "pairwise_table.csv", emit: pairwise_table_csv

    script:
    """
    pairwise_table_fasta.py variants_in_fasta.csv reformatted.fa
    #pairwise_table_fasta_downsample.py variants_in_fasta.csv reformatted.fa
    #pairwise_table_fasta_sng_end_restrict.py variants_in_fasta.csv reformatted.fa ${single_end_restrict}
    """
}


process PAIRWISE_BIALLELIC_TABLE {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path pairwise_table_csv
        val path_fn_modifier

    output:
        path "pairwise_biallelic_table.csv", emit: pairwise_biallelic_table_csv

    script:
    // NOTE
    // There was serious error which prevented it from working properly, but is now fixed.
    // Need to fix in other projects that use this also
    """
    biallelic_filter.py pairwise_table.csv
    """
}


process LOOKUP_TABLE_LDPOP {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val sample_size
        val theta
        val path_fn_modifier

    output:
        path "lookupTable.txt", emit: lookupTable_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores 4 -n ${sample_size} -th ${theta} -rh ${params.ldpop_rho_range} --approx > lookupTable.txt
    """
}


// process PYRHO_HAP_SETS_AND_MERGE {
//     publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

//     maxForks 1

//     input:
//         path lookup_table_txt
//         path pairwise_biallelic_table_csv
//         path seqgenOut
//         val genome_size
//         val path_fn_modifier

//     output:
//         path "table_ids_for_eq3.csv", emit: table_ids_for_eq3_csv
//         path "eq3.csv", emit: eq3_csv

//     script:
//     """
//     pyrho_hap_sets_and_merge.py lookupTable.txt pairwise_biallelic_table.csv seqgenOut.fa ${params.ldpop_rho_range} ${genome_size} > table_ids_for_eq3.csv
//     """
// }


process PAIRWISE_LOOKUP_FORMAT {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path pairwise_biallelic_table_csv
        val path_fn_modifier

    output:
        path "lookup_format.csv", emit: lookup_format_csv

    script:
    """
    pairwise_lookup_format_pyrho.py pairwise_biallelic_table.csv
    """
}


process CUSTOM_HAP_SETS_AND_MERGE {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path lookup_table_txt
        path pairwise_biallelic_table_csv
        path lookup_format_csv
        val ldpop_rho_range
        val path_fn_modifier

    output:
        path "table_ids_for_eq3.csv", emit: table_ids_for_eq3_csv
        path "eq3.csv", emit: eq3_csv

    script:
    """
    custom_hap_sets_and_merge.py lookupTable.txt pairwise_biallelic_table.csv lookup_format.csv ${ldpop_rho_range} > table_ids_for_eq3.csv
    """
}


process WATTERSON_ESTIMATE {
    maxForks 1

    input:
        val sample_size
        val genome_size
        path variants_in_fasta_csv

    output:
        stdout emit: theta

    script:
    """
    watterson_estimate.py variants_in_fasta.csv ${genome_size} ${sample_size}
    """
}


process P_IJ_GRID {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path eq3_csv
        val genome_size
        val recom_tract_len
        val ldpop_rho_range
        val path_fn_modifier

    output:
        path "p_ij_grid.csv", emit: p_ij_grid_csv

    script:
    """
    pij_grid_vectorised.py ${genome_size} ${recom_tract_len} ${ldpop_rho_range} eq3.csv
    """
}


process PAIRWISE_ESTIMATOR {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // errorStrategy 'ignore'

    // echo true
    
    maxForks 1
    
    input:
        path eq3_csv
        path p_ij_grid_csv
        path table_ids_for_eq3_csv
        path lookup_table
        val ldpop_rho_range
        val path_fn_modifier

    output:
        path "collected_likelihoods.csv", emit: collected_likelihoods_csv
    
    script:
    """
    pairwise_rho_estimator_intp_rect_biv.py eq3.csv table_ids_for_eq3.csv p_ij_grid.csv lookupTable.txt ${ldpop_rho_range}
    """

}


process FINAL_RESULTS {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // echo true

    maxForks 1

    input:
        path collectedFile
        val theta
        val path_fn_modifier

    output:
        path "final_results.txt", emit: final_results_txt

    script:
    """
    final_results.py collected_likelihoods.csv ${theta}
    """
}


process PROCESS_OUTPUT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path final_results_txt
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "processed_results.csv", emit: processed_results_csv

    script:
        """
        custom_est_process_output.py final_results.txt ${rho_rate} ${sample_size} ${genome_size}
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

    params.seed = 123
    params.mutation_rate = 0.01
    params.recom_tract_len = 500
    params.ldpop_rho_range = "101,100"
    params.effective_pop_size = 1
    params.rho_rates = 0.05
    params.sample_sizes  = 20
    params.genome_sizes = 25000

    params.single_end_restrict = 150
    
    // precomputed likelihood table
    // lookup_Table = Channel.fromPath("$baseDir/lookupTable.txt")
    
    // trees = Channel.fromPath("$baseDir/trees.txt")
    
    RATE_SELECTOR(params.rho_rates, params.mutation_rate, params.genome_sizes, params.sample_sizes, params.seed)

    // MS(RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, params.seed, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    FAST_SIM_BAC(RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.theta, RATE_SELECTOR.out.seed, RATE_SELECTOR.out.path_fn_modifier)

    // MS_PRIME(RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(MS.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(trees, RATE_SELECTOR.out.path_fn_modifier)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.seed, RATE_SELECTOR.out.path_fn_modifier)

    // ALLELE_FREQ(SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.path_fn_modifier)

    REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.path_fn_modifier)

    FASTA_VARIANT_SITES(REFORMAT_FASTA.out.reformatted_fa, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_TABLE(FASTA_VARIANT_SITES.out.variants_in_fasta_csv, REFORMAT_FASTA.out.reformatted_fa, params.single_end_restrict, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_BIALLELIC_TABLE(PAIRWISE_TABLE.out.pairwise_table_csv, RATE_SELECTOR.out.path_fn_modifier)

    LOOKUP_TABLE_LDPOP(RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.theta, RATE_SELECTOR.out.path_fn_modifier)

    // PYRHO_HAP_SETS_AND_MERGE(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv, SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_LOOKUP_FORMAT(PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv, RATE_SELECTOR.out.path_fn_modifier)

    CUSTOM_HAP_SETS_AND_MERGE(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv, PAIRWISE_LOOKUP_FORMAT.out.lookup_format_csv, params.ldpop_rho_range, RATE_SELECTOR.out.path_fn_modifier)

    WATTERSON_ESTIMATE(RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, FASTA_VARIANT_SITES.out.variants_in_fasta_csv)

    // P_IJ_GRID(PYRHO_HAP_SETS_AND_MERGE.out.eq3_csv, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    P_IJ_GRID(CUSTOM_HAP_SETS_AND_MERGE.out.eq3_csv, RATE_SELECTOR.out.genome_size, params.recom_tract_len, params.ldpop_rho_range, RATE_SELECTOR.out.path_fn_modifier)

    // PAIRWISE_ESTIMATOR(PYRHO_HAP_SETS_AND_MERGE.out.eq3_csv, PYRHO_HAP_SETS_AND_MERGE.out.table_ids_for_eq3_csv, P_IJ_GRID.out.p_ij_grid_csv, LOOKUP_TABLE_LDPOP.out.lookupTable_txt, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_ESTIMATOR(CUSTOM_HAP_SETS_AND_MERGE.out.eq3_csv, CUSTOM_HAP_SETS_AND_MERGE.out.table_ids_for_eq3_csv, P_IJ_GRID.out.p_ij_grid_csv, LOOKUP_TABLE_LDPOP.out.lookupTable_txt, params.ldpop_rho_range, RATE_SELECTOR.out.path_fn_modifier)

    FINAL_RESULTS(PAIRWISE_ESTIMATOR.out.collected_likelihoods_csv, WATTERSON_ESTIMATE.out.theta, RATE_SELECTOR.out.path_fn_modifier)

    PROCESS_OUTPUT(FINAL_RESULTS.out.final_results_txt, RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // collectedFile = PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    // PLOT_RESULTS(collectedFile)
}