# CDCgov/mycosnp-nf: Parameters

```bash
------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  CDCgov/mycosnp-nf v1.6.0
------------------------------------------------------
Typical pipeline command:

  nextflow run CDCgov/mycosnp-nf -profile singularity,test

Workflow selection options
  --workflow                   [string]  Name of the workflow to run (`PRE_MYCOSNP` for Pre-MycoSNP workflow | `NFCORE_MYCOSNP` for main MycoSNP workflow). 
                                         [default: NFCORE_MYCOSNP] 

Input/output options
  --input                      [string]  Path to comma-separated file containing information about the samples in the experiment.
  --add_sra_file               [string]  Path to comma-separated file containing SRA ids to download from NCBI. Format: Name,SRAID
  --add_vcf_file               [string]  Path to .csv (only one column) containing a list of file paths to vcf files generated from previous runs of this 
                                         workflow to include in this analysis. They must use the same exact reference. *.tbi file must be in same location. Each 
                                         line of .csv file has format: /path/to/vcf/file.gz 
  --outdir                     [string]  Path to the output directory where the results will be saved. [default: ./results]
  --publish_dir_mode           [string]  Method used to save pipeline results to output directory. [default: copy]

Reference genome options
  --fasta                      [string]  Path to FASTA formatted reference genome file.
  --ref_dir                    [string]  Path to reference genome masked files/picard dict/samtools fai/bwa index from previous run. If you use this command, it 
                                         invalidates `--fasta --ref_masked_fasta --ref_fai --ref_bwa --ref_dict` 
  --ref_masked_fasta           [string]  Path to reference genome masked fasta file. If you use this command, must provide all `--ref_masked_fasta --ref_fai 
                                         --ref_bwa --ref_dict` 
  --ref_fai                    [string]  Path to reference genome samtools fai file. If you use this command, must provide all `--ref_masked_fasta --ref_fai 
                                         --ref_bwa --ref_dict` 
  --ref_dict                   [string]  Path to reference genome picard tools dict file. If you use this command, must provide all `--ref_masked_fasta 
                                         --ref_fai --ref_bwa --ref_dict` 
  --ref_bwa                    [string]  Path to reference genome bwatools bwa directory. If you use this command, must provide all `--ref_masked_fasta 
                                         --ref_fai --ref_bwa --ref_dict` 
  --snpeffconfig               [string]  Path to snpeff config [default: ${projectDir}/assets/snpeffdb/]
  --species                    [string]  Species name [default: candida_auris_gca_016772135.1]
  --positions                  [string]  FKS1 hotspot coordinates. Non-synonymous variants in these regions will be included in the snpeffr report 
                                         (snpeff/combined_cauris_refB11205_fks1.csv). [default: 
                                         fks1_hs1=221637:221663,fks1_hs2=223782:223805,fks1_hs3=221805:221807] 
  --genes                      [string]  Genes [default: CAB11_002014]
  --exclude                    [string]  exclude variants [default: synonymous_variant]

MycoSNP Run/Save/Skip Options
  --save_reference             [boolean] Saves the reference genome/index files to the results folder [default: true]
  --save_alignment             [boolean] Saves the reference alignment BAM files to the samples results folder [default: true]
  --rapidnj                    [boolean] Build a tree using the RapidNJ neighbour-joining algorithm [default: true]
  --fasttree                   [boolean] Build a tree using the FastTree approximate ML algorithm [default: true]
  --snpeff                     [boolean] Run snpEff on the final VCF file [default: false]
  --iqtree                     [boolean] Build a tree using the IQ-TREE ML algorithm
  --raxmlng                    [boolean] Build a tree using the RAxML-NG ML algorithm
  --save_debug                 [boolean] Save intermediate files for debugging
  --skip_samples               [string]  a comma separated list of IDs to skip for variant calling and analysis
  --skip_samples_file          [string]  a file with new-line separated list of IDs to skip for variant calling and analysis
  --skip_combined_analysis     [boolean] Skip combined variant analysis (run reference prep and mapping)
  --skip_phylogeny             [boolean] Skip phylogenetic tree creation

MycoSNP Run Params
  --sample_ploidy              [integer] Ploidy of sample (GATK) [default: 1]
  --coverage                   [integer] Coverage is used to calculate a down-sampling rate that results in the specified coverage. For example, if coverage is 
                                         70, then FASTQ files are down-sampled such that, when aligned to the reference, the result is approximately 70x 
                                         coverage [default: 0] 
  --gvcfs_filter               [string]  Filter criteria for variants (GATK) [default: QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10]
  --gatkgenotypes_filter       [string]  Filter criteria for script filterGatkGenotypes [default: --min_GQ "50" --keep_GQ_0_refs --min_percent_alt_in_AD 
                                         "0.8" --min_total_DP "10" --keep_all_ref] 
  --max_amb_samples            [integer] Max number of samples with ambiguous calls for inclusion (GATK) [default: 10000000]
  --max_perc_amb_samples       [integer] Max percent of samples with ambiguous calls for inclusion (GATK) [default: 10]
  --min_depth                  [integer] Min depth for a base to be called as the consensus sequence, otherwise it will be called as an N. Set to 0 to 
                                         disable. [default: 10] 
  --mask                       [boolean] Perform masking of reference genome before analysis [default: true]

Pre-MycoSNP options
  --assembler                  [string]  Assembler to use with Shovill. Options are 'skesa', 'velvet', 'megahit', 'spades' [default: skesa]
  --shovill_depth              [integer] Downsample FASTQ files to this depth before assembly [default: 70]
  --genome_size                [string]  The approx. genome size of the sample, used by Shovill. Will be automatically determined if left blank.
  --min_contig_cov             [integer] Minimum contig depth of coverage for it to be retained in the assembly, used by Shovill [default: 10]
  --min_contig_len             [integer] Minimum contig length for it to be retained in the assembly, used by Shovill [default: 300]
  --gambit_h5                  [string]  Path to the GAMBIT '.h5' file. [default: 
                                         gs://theiagen-public-files-rp/terra/theiaeuk-files/gambit/221130-theiagen-fungal-v0.2.h5] 
  --gambit_db                  [string]  Path to the GAMBIT '.db' file. [default: 
                                         gs://theiagen-public-files-rp/terra/theiaeuk-files/gambit/221130-theiagen-fungal-v0.2.db] 
  --subtype_db                 [string]  Path to a directory containing the subtyper files. This directory should contain sourmash signature files, each 
                                         containing multiple sketches for the representative subtypes. This directory should also contain a csv file called 
                                         sourmash_taxa.csv mapping each taxon name to a sourmash signature file. Taxon names must be the same as what is 
                                         reported by GAMBIT. [default: ${projectDir}/assets/sourmash_db/] 

Generic options
  --help                       [boolean] Display help text.
  --tmpdir                     [string]  temporary directory for shell and java processes [default: /tmp]
```