

## Create Conda ENV with Nextflow & NF-Core
```
mamba env create -n nextflow -c bioconda -c conda-forge nf-core nextflow git graphviz yamllint
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
nf-core modules install rapidnj
nf-core modules install raxmlng
nf-core modules install samtools/view
```
