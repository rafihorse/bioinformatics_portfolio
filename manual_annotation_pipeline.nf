log.info """\n
    PACBIO MANUAL ANNOTATION NEXTFLOW PIPELINE
    =====================================================
    genome: ${params.genome}
    reads: ${params.fastq}
    locus: ${params.locus}
    gtf: ${params.gtf}
    biomart: ${params.biomart}
    outdir: ${params.outdir}
    runPart2: ${params.runPart2}
    big: ${params.big}
    """
    .stripIndent(true)

// Check if 'chr' is present in params.locus
def hasChrPrefix = params.locus.contains('chr')
log.info "Locus contains 'chr' prefix: ${hasChrPrefix}"

process INDEX {
    label 'small'

    input:
    path genome

    output:
    path "${genome}.fai"

    script:
    """
    samtools faidx $genome
    """
}

process CHUNK {
    label 'medium'
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("${fastq}.split/${fastq.simpleName}.part_*.fastq.gz")

    script:
    """
    seqkit split2 --by-size 1000000 -e .gz -j ${task.cpus} $fastq
    """
}

process MAP {
    label 'medium'
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    path genome
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("${fastq.baseName.replaceAll('.fastq', '')}.sam")

    script:
    """
    minimap2 -ax splice:hq -I 8g --MD -t ${task.cpus} $genome $fastq -o ${fastq.baseName.replaceAll('.fastq', '')}.sam
    """
}

process SORT {
    label 'small'

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path("${sam.baseName}.bam")

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sam.baseName}.bam -O BAM -m 3G $sam
    """
}

process FLAGSTAT {
    label 'small'
    publishDir "stats/flagstat", mode: 'copy'

    input:
    tuple val(id), path(bam)

    output:
    path "${bam.simpleName}_flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${bam.simpleName}_flagstat.txt
    """
}

process INDEX_BAM {
    label 'medium'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path(bam), path("${bam}.bai")

    script:
    """
    samtools index -b -@ ${task.cpus} ${bam} --output ${bam}.bai
    """
}

process CALL_VCF {
    label 'medium'

    input:
    path(bams)
    path(index)
    path(reference)
    path(fai)

    output:
    tuple path("all_samples_minfrac0.2_${params.locus.replaceAll(",","")}.vcf.gz"), path("all_samples_minfrac0.2_${params.locus.replaceAll(",","")}.vcf.gz.tbi")

    script:
    """
    freebayes --bam ${bams} --region ${params.locus} --fasta-reference ${reference} \
    --pooled-discrete --min-coverage 10 \
    --min-alternate-count 10 \
    --min-alternate-fraction 0.1 \
    --max-complex-gap 40 \
    --haplotype-length 0 > all_samples_minfrac0.2_${params.locus.replaceAll(",","")}.vcf

    bgzip all_samples_minfrac0.2_${params.locus.replaceAll(",","")}.vcf
    tabix -p vcf all_samples_minfrac0.2_${params.locus.replaceAll(",","")}.vcf.gz
    """

}

process CORRECT_BAM {
    label 'medium'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/seqs", mode: 'copy'
    
    input:
    tuple val(id), path(bam), path(bai)
    path(genome)
    tuple path(vcf), path(vcf_index)

    output:
    tuple val(id), path("${id}_clean.bam"), emit: bam
    tuple val(id), path("${id}_clean.bam.bai"), emit: bai
    tuple val(id), path("${id}_clean.fa"), emit: fasta

    script:
    """
    samtools view -@ ${task.cpus} -h -o ${id}.sam ${bam} #${params.locus}
    transcriptclean -s ${id}.sam -g ${genome} -t ${task.cpus} -v ${vcf} --maxLenIndel 5 --correctSJs=false --primaryOnly --deleteTmp -o ${id}
    samtools view -@ ${task.cpus} -b -o ${id}_clean.bam ${id}_clean.sam
    samtools index -b -@ ${task.cpus} ${id}_clean.bam --output ${id}_clean.bam.bai
    """
}

process ANNOTATE_TRANSCRIPTS {
    label 'medium'
    publishDir "${params.outdir}/isoforms", mode: 'copy'
    cache 'false'
    
    input:
    path(bams)
    path(bais)
    path(reference)
    tuple path(vcf), path(vcf_index)

    output:
    path("all_bams_SJ_*-supporting_reads.bam"), emit: sj_bams
    path("all_bams_sj.bed12")
    path("all_bams_sj_counts.tsv")
    path("all_bams_sj_ucsc.bed12")
    path("all_bams_sj.fasta")
    path("all_bams_abundant_sj_counts.tsv")
    path("all_bams_abundant_sj.bed12"), emit: bed
    path("all_bams_abundant_sj_ucsc.bed12"), emit: ucsc
    path("all_bams_abundant_sj.fasta"), emit: fasta
    path("*_read_frequencies.png")
    path("*_start_positions.bw")
    path("*_stop_positions.bw")

    script:
    """
    manual_isoform_annotation_updated.py ${bams} --locus ${params.locus} --reference $reference --output all_bams
    """
    
}

process GET_FASTA {
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(genome)
    path(bed)

    output:
    path("manual_annotations.fa")

    script:
    """
    bedtools getfasta -fi ${genome} -bed ${bed} -name -s -split -fo manual_annotations.fa
    """
}

process GET_ORFS {
    label 'small'
    publishDir "${params.outdir}/orfs", mode: 'copy'

    input:
    path(fasta)

    output:
    path("manual_annotations.orfs.fa")

    script:
    """
    esl-translate -l 100 -m ${fasta} > manual_annotations.orfs.fa
    """
}

process ANNOTATE_CDS {
    label 'small'
    publishDir "${params.outdir}/orfs", mode: 'copy'

    input:
    path(bed)
    path(fasta)

    output:
    path("manual_annotations_cds.bed"), emit: bed
    path("manual_annotations_cds_ucsc.bed"), emit: ucsc

    script:
    """
    add_orfs_to_bed.py --bed_file ${bed} --fasta_file ${fasta} --output manual_annotations_cds.bed --min_orf_size 200
    """
}

process CREATE_GTF {
    label 'small'
    publishDir "${params.outdir}/transcriptome", mode: 'copy'

    input:
    path(bed)
    path(biomart)

    output:
    path("manual_annotations_corrected.gtf")

    script:
    """
    gffread -T ${bed} -o manual_annotations_temp.gtf
    correct_gtf_attributes.py -g manual_annotations_temp.gtf -m ${biomart} -o manual_annotations_corrected.gtf -t ${launchDir}/manual_annotation/orfs/transcript_types.txt
    """
}

process CREATE_TRANSCRIPTOME {
    label 'small'
    publishDir "${params.outdir}/transcriptome", mode: 'copy'

    input:
    path(sra_gtf)
    path(gtf)
    path(genome)
    path(sra_fasta)

    output:
    path("full_manual_annotations_sorted.gtf")
    path("manual_annotations_transcriptome.fa")

    script:
    """
    # Remove SRA1 from the original GTF
    awk -F '\\t' 'BEGIN {IGNORECASE=1} \$9 !~ /gene_name "SRA1"/' ${gtf} > "${gtf.simpleName}_SRA1_removed.gtf"

    # Merge
    cat "${gtf.simpleName}_SRA1_removed.gtf" ${sra_gtf} > full_manual_annotations.gtf

    # Sort
    gffread full_manual_annotations.gtf -E -F -T -o full_manual_annotations_sorted.gtf
    gffread -w SRA1_removed_transcriptome.fa -g ${genome} ${gtf.simpleName}_SRA1_removed.gtf
    cat SRA1_removed_transcriptome.fa ${sra_fasta} > manual_annotations_transcriptome.fa
    """
}

process GENOME_BROWSER {
    label 'small'
    publishDir "${params.outdir}/orfs", mode: 'copy'

    input:
    path(bed)

    output:
    path("manual_annotations_orfs_browser.png")

    script:
    """
    genome_browser.py -i ${bed} -o manual_annotations_orfs_browser.png
    """
}

process MERGE_BAM {
    label 'medium'

    input:
    tuple val(id), path(bam_files)

    output:
    tuple val(id), path("${id}_merged.bam")

    script:
    """
    samtools merge -@ ${task.cpus} -o ${id}_merged.bam ${bam_files.join(' ')}
    """

}

workflow {
    main:
        Channel
            .fromPath(params.genome)
            .set { genome_file }

        // FIX: handle either (a) glob of FASTQ files or (b) a single list file
        def fastq_pattern = params.fastq
        def fastq_ch = Channel.fromPath(fastq_pattern)
                              .flatMap { f ->
                                  // If it's a single list file of paths (txt/list), expand; else pass file itself
                                  if (f.name =~ /(?i)\.(txt|list)$/ && f.size() < 5_000_000) {
                                      f.text.readLines()
                                       .findAll { it && !it.startsWith('#') }
                                       .collect { file(it) }
                                  } else {
                                      [f]
                                  }
                              }

        fastq_ch
            .map { f -> [f.simpleName, f] }
            .groupTuple()
            .set { fastq_files }

        Channel
            .fromPath(params.gtf)
            .set { gtf_file }

        Channel
            .fromPath(params.biomart)
            .set { biomart_file }

        index_ch = INDEX(genome_file)

        // fastq_files = fastq_files.randomSample(100, 42)

        // New workflow branch: split & map chunks if "big", otherwise map fastq files directly
        if (params.big) {
            chunk_ch   = CHUNK(fastq_files)
            // Flatten the chunk channel from: [id, [chunk1, chunk2, ...]]
            flat_chunks = chunk_ch
                .map { id, chunks -> [id, chunks instanceof List ? chunks : [chunks]] }
                .flatMap { id, chunks -> chunks.collect { [id, it] } }
            //flat_chunks.view()
            map_chunks = MAP(genome_file.first(), flat_chunks)
            sort_chunks = SORT(map_chunks)
            // Group sorted chunks by id and extract bam files
            grouped_chunks = sort_chunks.groupTuple()
            merged = MERGE_BAM(grouped_chunks)
            sort_ch = merged
        } else {
            map_ch  = MAP(genome_file.first(), fastq_files)
            sort_ch = SORT(map_ch)
        }

        flagstat_ch   = FLAGSTAT(sort_ch)
        index_bam_ch  = INDEX_BAM(sort_ch)

        vcf_ch = CALL_VCF(index_bam_ch.map { it[1] }.collect(), index_bam_ch.map { it[2] }.collect(), genome_file, index_ch)
        bam_ch = CORRECT_BAM(index_bam_ch, genome_file.first(), vcf_ch.first())

        // Extract just the BAM files and BAI files without the sample IDs
        bam_files_only = bam_ch.bam.map { id, file -> file }
        bai_files_only = bam_ch.bai.map { id, file -> file }

        annotated_ch = ANNOTATE_TRANSCRIPTS(
            bam_files_only.collect(),
            bai_files_only.collect(),
            genome_file.first(),
            vcf_ch.first()
        )

        orfs_ch = GET_ORFS(annotated_ch.fasta)

        cds_bed = ANNOTATE_CDS(annotated_ch.bed, orfs_ch)
        
        // Use UCSC format bed if locus contains 'chr', otherwise use standard bed
        if (hasChrPrefix) {
            browser = GENOME_BROWSER(cds_bed.ucsc)
        } else {
            browser = GENOME_BROWSER(cds_bed.bed)
        }

        if (params.runPart2) {
            corrected_gtf = hasChrPrefix ? CREATE_GTF(cds_bed.ucsc, biomart_file) : CREATE_GTF(cds_bed.bed, biomart_file)
            transcriptome_ch = CREATE_TRANSCRIPTOME(corrected_gtf, gtf_file, genome_file, annotated_ch.fasta)
        }
}