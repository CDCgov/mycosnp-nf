

## Create Conda ENV with Nextflow & NF-Core
```
mamba create -n nextflow -c bioconda -c conda-forge nf-core nextflow git graphviz yamllint
conda activate nextflow
```

## NF-Core Create New pipeline from template
```
nf-core create
```

>```
>Workflowname: mycosnp
>Description: From Mycosnp
>Author: CDC
>```

## Add GitLab remote repo
```
cd nf-core-mycosnp
git remote add origin git@git.biotech.cdc.gov:scbs/nextflow/mycosnp.git
git push --all origin
```

## Create new nf-core modules

```
conda activate nextflow
```

Before you start, please check that the module you wish to add isn't already on [nf-core/modules](https://github.com/nf-core/modules.git):

* Use the `nf-core modules list` command
* Check [open pull requests](https://github.com/nf-core/modules/pulls)
* Search [open issues](https://github.com/nf-core/modules/issues)

If the module doesn't exist on nf-core/modules:

* Please create a [new issue](https://github.com/nf-core/modules/issues/new?assignees=&labels=new%20module&template=new_nodule.md&title=new%20module:) before adding it
* Set an appropriate subject for the issue e.g. `new module: fastqc`
* Add yourself to the `Assignees` so we can track who is working on the module

* Check that the new module you've added follows the [new module guidelines](https://nf-co.re/developers/adding_modules#new-module-guidelines-and-pr-review-checklist)

* Lint the module locally to check that it adheres to nf-core guidelines before submission

```bash

```


```bash
nf-core modules install bwa/index
nf-core modules install bwa/mem
nf-core modules install mummer
nf-core modules install nucmer
nf-core modules install bedtools/maskfasta
nf-core modules install samtools/faidx
nf-core modules install samtools/index
nf-core modules install samtools/sort
nf-core modules install gatk4/haplotypecaller
nf-core modules install gatk4/genotypegvcfs
nf-core modules install gatk4/variantfiltration
nf-core modules install bcftools/index
nf-core modules install bcftools/query
nf-core modules install bcftools/consensus
nf-core modules install fastqc
nf-core modules install multiqc
nf-core modules install qualimap/bamqc
nf-core modules install picard/markduplicates
nf-core modules install seqtk/sample
nf-core modules install bcftools/view
nf-core modules install gatk4/selectvariants
nf-core modules install gatk4/combinegvcfs
nf-core modules install picard/createsequencedictionary
nf-core modules install samtools/view
nf-core modules install picard/cleansam
nf-core modules install picard/fixmateinformation
nf-core modules install picard/addorreplacereadgroups
nf-core modules install seqtk/seq
nf-core modules install seqtk/rename
nf-core modules install faqcs
nf-core modules install seqkit/pair
nf-core modules install seqkit/replace
nf-core modules install rapidnj
nf-core modules install raxmlng
nf-core modules install fasttree
nf-core modules install iqtree
nf-core modules install samtools/view
nf-core modules install samtools/stats
nf-core modules install samtools/idxstats
nf-core modules install samtools/flagstat
nf-core modules install snpdists
```

```bash
nf-core modules list local
```

```bash
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.2



INFO     Modules installed in '.':                                                                                                                     list.py:128
                                                                                                                                                                  
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
┃ Module Name                     ┃ Repository      ┃ Version SHA                              ┃ Message                                            ┃ Date       ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
│ bcftools/consensus              │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ bcftools/index                  │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ bcftools/query                  │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ bcftools/view                   │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ bedtools/maskfasta              │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ bwa/index                       │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ bwa/mem                         │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ custom/dumpsoftwareversions     │ nf-core/modules │ 20d8250d9f39ddb05dfb437603aaf99b5c0b2b41 │ Update all modules to new NF DSL2 syntax (#1099)   │ 2021-11-26 │
│ faqcs                           │ nf-core/modules │ 2cd502a236aec1a1eecbe1f9189af3414efbf1d8 │ Faqcs patch (#1367)                                │ 2022-03-02 │
│ fastqc                          │ nf-core/modules │ 9d0cad583b9a71a6509b754fdf589cbfbed08961 │ Change syntax from task.ext.suffix to              │ 2021-12-02 │
│                                 │                 │                                          │ tast.ext.prefix in all modules (#1110)             │            │
│ fasttree                        │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ gatk4/combinegvcfs              │ nf-core/modules │ a25423dbb974e2e69268b5cdbc12e5bc1558e053 │ Gatk4 combinegvcfs (#1342)                         │ 2022-02-23 │
│ gatk4/genotypegvcfs             │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ gatk4/haplotypecaller           │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ gatk4/selectvariants            │ nf-core/modules │ 640031762334fd1d534d667c2e1a982c5513a7aa │ Gatk4 selectvariants (#1346)                       │ 2022-02-24 │
│ gatk4/variantfiltration         │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ iqtree                          │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ multiqc                         │ nf-core/modules │ 20d8250d9f39ddb05dfb437603aaf99b5c0b2b41 │ Update all modules to new NF DSL2 syntax (#1099)   │ 2021-11-26 │
│ mummer                          │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ nucmer                          │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ picard/addorreplacereadgroups   │ nf-core/modules │ a0d91e4a93daab335e8dbae31b345188c120e0a0 │ Picard addorreplacereadgroups (#1305)              │ 2022-02-19 │
│ picard/cleansam                 │ nf-core/modules │ 8a20253f4028133b589c2c169aaa68a5a7fe848d │ update args & convert to bam (#1355)               │ 2022-02-25 │
│ picard/createsequencedictionary │ nf-core/modules │ 62e5d1f0b38367d628e6b4c32b45d462f2c0b325 │ Picard createsequencedictionary (#1310)            │ 2022-02-19 │
│ picard/fixmateinformation       │ nf-core/modules │ f655e5dea2e25f403e23fb6dfdc33cc538666769 │ Picard fixmateinformation (#1315)                  │ 2022-02-19 │
│ picard/markduplicates           │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ qualimap/bamqc                  │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ rapidnj                         │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ raxmlng                         │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ samtools/faidx                  │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ samtools/flagstat               │ nf-core/modules │ 1ad73f1b2abdea9398680d6d20014838135c9a35 │ update samtools version to 1.15 (#1358)            │ 2022-02-28 │
│ samtools/idxstats               │ nf-core/modules │ 1ad73f1b2abdea9398680d6d20014838135c9a35 │ update samtools version to 1.15 (#1358)            │ 2022-02-28 │
│ samtools/index                  │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ samtools/sort                   │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ samtools/stats                  │ nf-core/modules │ 1ad73f1b2abdea9398680d6d20014838135c9a35 │ update samtools version to 1.15 (#1358)            │ 2022-02-28 │
│ samtools/view                   │ nf-core/modules │ e745e167c1020928ef20ea1397b6b4d230681b4d │ Fix formatting in yaml files, add yamllint config  │ 2022-02-15 │
│                                 │                 │                                          │ (#1279)                                            │            │
│ seqkit/pair                     │ nf-core/modules │ 4c59984d7bcc772598ba44fb8289051352c7a6c6 │ Seqkit pair (#1348)                                │ 2022-02-24 │
│ seqkit/replace                  │ nf-core/modules │ 24f0bdd14ec32e0114aa6ee5337ddbd490ffd570 │ added module seqkit replace (#1382)                │ 2022-03-09 │
│ seqtk/rename                    │ nf-core/modules │ a69faefee8ff63565629a2b6dfc40e08f4690b80 │ Seqtk rename (#1304)                               │ 2022-02-16 │
│ seqtk/sample                    │ nf-core/modules │ e20e57f90b6787ac9a010a980cf6ea98bd990046 │ Add when: block (#1261)                            │ 2022-02-04 │
│ seqtk/seq                       │ nf-core/modules │ 1016c9bd1a7b23e31b8215b8c33095e733080735 │ Seqtk seq (#1340)                                  │ 2022-02-23 │
└─────────────────────────────────┴─────────────────┴──────────────────────────────────────────┴────────────────────────────────────────────────────┴────────────┘
```


## Bumping a pipeline version number

When releasing a new version of a nf-core pipeline, version numbers have to be updated in several different places. The helper command nf-core bump-version automates this for you to avoid manual errors (and frustration!)

```
nf-core bump-version 1.2
```