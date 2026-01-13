process VCF_QC {
    tag "vcf-qc"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    path(vcffasta)

    output:
    path("vcf-qc-report.txt"),        emit: vcf_qc_report
    path("vcfqc_passed_samples.txt"), emit: vcf_qc_passed
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Prepare main QC report file
    printf "Sample Name\\tLength\\tNumber-N\\t%%N\\n" > vcf-qc-report.txt

    # Process FASTA and generate reports
    awk -v threshold=$params.percent_n '
    BEGIN { passed_samples_count = 0 }

    /^>/ {
        if (NR > 1) {
            # Check for division by zero
            if (c > 0) {
                percent_n = (d * 100 / c)
            } else {
                percent_n = "NA"
            }
            printf "%s\\t%s\\t%.2f\\n", c, d, percent_n >> "vcf-qc-report.txt"

            if (percent_n != "NA" && percent_n < threshold) {
                print sample_name >> "vcfqc_passed_samples_ref.txt"
                if (passed_samples_count > 0) {
                    print sample_name >> "vcfqc_passed_samples.txt"
                }
            }
            passed_samples_count++
        }
        c = 0; d = 0; sample_name = substr(\$0, 2)
        printf "%s\\t", sample_name >> "vcf-qc-report.txt"
    }

    /^[^>]/ {
        c += length(\$0)
        d += gsub(/[Nn]/, "")
    }

    END {
        # Process the last sample
        if (c > 0) {
            percent_n = (d * 100 / c)
        } else {
            percent_n = "NA"
        }
        printf "%s\\t%s\\t%.2f\\n", c, d, percent_n >> "vcf-qc-report.txt"

        if (percent_n != "NA" && percent_n < threshold) {
            print sample_name >> "vcfqc_passed_samples_ref.txt"
            if (passed_samples_count > 0) {
                print sample_name >> "vcfqc_passed_samples.txt"
            }
        }

        # If no samples passed, ambiguity threshold was exceeded. Log a warning
        if (passed_samples_count == 0) {
            printf "All SNPs were filtered because the ambiguity threshold was exceeded.\\n" >> "vcf-qc-report.txt"
        }
    }' ${vcffasta}

    # Ensure all outputs are generated
    if [ ! -s "vcfqc_passed_samples.txt" ]; then
        touch vcfqc_passed_samples.txt
    fi

    if [ ! -s "vcfqc_passed_samples_ref.txt" ]; then
        touch vcfqc_passed_samples_ref.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/^GNU Awk //; s/^mawk //; s/,.*//; s/ .*//; s/awk version //')
    END_VERSIONS
    """

    stub:
    """
    touch vcf-qc-report.txt
    touch vcfqc_passed_samples.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/^GNU Awk //; s/^mawk //; s/,.*//; s/ .*//; s/awk version //')
    END_VERSIONS
    """
}
