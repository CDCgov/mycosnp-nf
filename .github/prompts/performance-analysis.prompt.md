---
agent: agent
model: Claude Sonnet 4
---

## Context and Goal

The goal of this analysis is to evaluate the performance of a bioinformatics pipeline for performing
whole genome sequencing analysis of fungal organisms (e.g. _Candida auris_) from Illumina paired-end reads,
specifically focusing on its efficiency, resource utilization, and identification of bottlenecks.
The pipeline involves several key steps:

- Sequencing reads quality metrics
- Remove unpaired reads
- Trim/filter reads
- De novo assembly (downsampling and assembly)
- Taxonomic classification (classifies isolate to genus/species level)
- Subtyping (by default, for _C. auris_, Pre-MycoSNP performs clade typing)
- Summary report containing genus/species classification, subtype, and read and assembly quality metrics

## Input Data

Attached or included execution trace file.

## Performance Metrics

The performance analysis should consider the following key metrics:

- **Execution Time:** Total time taken for each stage of the pipeline and the overall workflow.
- **Resource Utilization:** CPU usage, memory consumption, and disk I/O at each stage.
- **Error Rate:** Incidence of errors or failures within the pipeline steps.

## Analysis Request

The performance report should be presented in a clear and structured format, potentially using markdown for easier readability and stored in the 'docs' folder. It should include:

- **Identify Bottlenecks:** Pinpoint the specific steps or tasks within the pipeline that are consuming the most resources or causing delays.
- **Suggest Optimization Strategies:** Provide recommendations for improving the pipeline's performance, including potential changes to parameters, algorithms, or infrastructure.
- **Generate Performance Report:** Create a summary report highlighting the key findings, including visualizations (e.g., graphs of resource utilization over time, bottleneck identification) if possible.

## Example Output Format (Desired)

The performance report should be presented in a clear and structured format, potentially using markdown for easier readability and stored in the 'docs' folder. It should include:

- An executive summary outlining the main findings.
- Detailed sections for each performance metric, including numerical data and explanations.
- Recommendations for optimization, prioritized by potential impact.
- Provide targeted code or config level suggestions.
- Provide targeted resource (CPUs and memory) adjustments or optimizations based on observed usage. Be conservative in cpu and memory adjustments
