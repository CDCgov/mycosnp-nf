---
applyTo: "**"
excludeAgent: ["coding-agent"]
---

# ph-core Code Review Standards

This file contains instructions for GitHub Copilot code review to ensure AMDP prerequisites and ph-core NextFlow standards compliance for bioinformatics workflow development.

## Review Scope

Apply these standards when reviewing:
- NextFlow workflow files (`.nf`)
- Module and subworkflow implementations
- Configuration files (`nextflow.config`, `*.config`)
- Documentation (README, CHANGELOG, meta.yml)
- Test files (nf-test)
- GitHub Actions workflows

## Compliance Tiers

### Tier 1: AMDP Prerequisites (MANDATORY)
Hard requirements for AMDP onboarding. Always check these first.

**Workflow Requirements:**
- ✅ Written in NextFlow with `.nf` files
- ✅ Single command execution (`nextflow run`) - pipelines must run with a single command
- ✅ Cloud compatibility (AWS/Azure/GCP executors) - must be tested on cloud computing environments
- ✅ Static assets and large test data stored separately (not in repo) and properly documented
- ✅ No deprecated NextFlow features (e.g., entry points replaced with profiles)
- ✅ Proper credits and acknowledgements
- ✅ Workflow specificity: single pipeline per data/analysis type

**Module Requirements:**
- ✅ Containers tagged in AMDP Elastic Container Registry
- ✅ Specific version tags (never `latest`)
- ✅ Pass Trivy security scans (no HIGH/CRITICAL vulnerabilities with no unresolved issues)

### Tier 2: ph-core Requirements
Additional standards for ph-core namespace workflows.

**Workflow:**
- ✅ Built with nf-core template (custom branded)
- ✅ MIT license
- ✅ Semantic versioning (e.g., v1.2.3)
- ✅ Git branches: `main`/`master`, `dev`, `TEMPLATE`
- ✅ Directory structure: `assets/`, `bin/`, `conf/`, `docs/`, `libs/`
- ✅ GitHub Actions CI/CD
- ✅ Passes `nf-core lint` tests
- ✅ nf-test coverage with minimal examples
- ✅ ph-core namespace identified in README
- ✅ SHARE IT Act metadata block in README
- ✅ Minimum inputs: pipelines should run with minimal required input
- ✅ Standardized parameters: strive for standardized usage patterns
- ✅ Docker support: software bundled using Docker with versioning
- ✅ Publication credit: acknowledge nf-core community in publications
- ✅ nf-core template version: created with nf-core/tools v3.2.0 or later
- ✅ Remove nf-core branding: select "Custom" during pipeline creation

**Modules:**
- ✅ Specific container versions (NOT `latest`)
- ✅ Single process per module
- ✅ Standard meta fields (`meta.id`, `meta.single_end`)
- ✅ Emit `versions.yml` (remove leading "v")
- ✅ Stored in `modules/` directory
- ✅ Include `meta.yml` descriptor file
- ✅ Include nf-test coverage

**Scripts:**
- ✅ Placed in `bin/` directory
- ✅ Shebang and executable permissions
- ✅ Called directly in processes

### Tier 3: ph-core Recommendations
Best practices for optimal quality.

- Modern file formats (CRAM over BAM)
- DOI assignment for citability
- Comprehensive test coverage
- Detailed documentation
- Code alignment standards ('Harshil Alignment™️')
- Optimized workflow size

**Subworkflows:**
- Minimum 2 modules per subworkflow
- Comments describing input/output channel structures
- Channel structure descriptions in meta.yml
- Test coverage with tags for dependent modules
- Use `assertAll()` with minimum success assertion and versions
- All output channels in nf-test snapshots or existence verification
- Test and assertions for each input/output type
- Descriptive test names (dataset and configuration)
- Named file extensions for outputs (exceptions for datasets/directories)
- Handle optional inputs with empty list `[]` workaround
- Lowercase directory naming
- `snake_case` for parameters and channels
- `camelCase` for functions
- Input channel names signify object type
- Output channel names based on major file type only
- No assumption of named params from parent workflow

## Key Standards Reference

### Naming Conventions
- **Workflows**: lowercase, no punctuation (`myworkflow`)
- **Modules**: lowercase directory (`modules/bwa/mem/`)
- **Processes**: UPPERCASE with underscore (`BWA_MEM`)
- **Parameters**: `snake_case`
- **Functions**: `camelCase`
- **Channels**: `snake_case`, lowercase
- **Outputs**: `${prefix}.ext` (e.g., `${prefix}.fq.gz`)

### Container Management
- Use specific versions: `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`
- Never use `latest` tag
- AMDP containers must be in Elastic Container Registry
- Run Trivy scans: `trivy image <container>`

### Module Structure Example
```nextflow
process TOOL_NAME {
    tag "$meta.id"
    label 'process_low'

    container 'container:specific_version'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tool $args \\
        -o ${prefix}.bam \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch versions.yml
    """
}
```

### Module Development Standards

**General Requirements:**
- All mandatory input files included in input channel definitions
- All optional input files included in input channel definitions
- Non-mandatory command-line tool non-file arguments provided via `$task.ext.args`
- Each command must have an `$args` variable; multiple args variables (`$args`, `$args2`, `$args3`) must correspond to order of tools in command line
- Process definition must not change the `when` statement
- Use shell `tee` command for capturing STDOUT and STDERR if needed
- Stub block must exist for all modules
- Module definitions with process blocks in order, no repeated `script:` blocks
- Non-file mandatory commands should be a value channel
- Limit multi-command piping unless there's runtime/storage advantage
- Inputs and outputs should be compressed where applicable
- If tool doesn't exit with code 0 on success, use `||` operator for alternative command (e.g., testing file size)
- All code aligned following 'Harshil Alignment™️' format
- Use only one tool per module if possible

**Input/Output Management:**
- Input channel path declarations for all possible input files (required and optional)
- Options/flags/parameters not required by tool should use `ext.args`, not `val` channel input
- Named file extensions for output channels (e.g., `path "*.txt", emit:txt`), with exceptions for datasets/directories
- Generic auxiliary files (e.g., reference sequences) may use dedicated input channel
- Input channel `val` declarations for all mandatory non-file inputs essential for tool function
- Optional inputs: use empty list `[]` workaround (not natively supported)
- Optional outputs marked as optional
- Each output file type emitted in own channel (max one), with meta map if provided (except `versions.yml`)

**Parameters:**
- All params within module only initialized and used in local context
- Parameters evaluated per sample (e.g., single-end/paired-end) defined within process
- Module file should only define input/output files as command-line parameters

**Resource Requirements:**
- Appropriate resource label: `process_single`, `process_low`, `process_medium`, or `process_high`
- If tool supports multi-threading: provide parameter using `$task.cpus`
- If no multi-threading: consider `process_single` unless large RAM required
- For multiple multi-threaded tools: assign CPUs per tool

**Software/Container Requirements:**
- Software requirements declared using Nextflow `container` directive
- Multi-tool containers available on BioContainers (e.g., bwa and samtools)
- New multi-tool containers: submit PR to BioContainers multi-package-containers repository

**Deprecation:**
- When alternative available on nf-core modules: add deprecation message at top and assert in code body

### Testing Requirements
- Use nf-test framework
- Minimal test data (subset datasets)
- Test all input/output types
- Include stub tests
- Use `assertAll()` for assertions
- One snapshot per test
- Descriptive test names

**Module-Specific Testing:**
- Tags for dependent modules must be specified to trigger tests when upstream modules change
- Minimum one success assertion and versions in snapshot
- Only one snapshot per module test containing all assertions
- Test each input/output type combination
- Test names describe dataset and configuration used
- Input data referenced with `modules_testdata_base_path` parameter

**Subworkflow-Specific Testing:**
- Tags for dependent modules ensure upstream changes trigger tests
- All output channels in snapshots or verified for existence
- Test each input/output type combination
- Descriptive test names indicating dataset and configuration

### Documentation Requirements

**README.md must include:**
- SHARE IT Act metadata block at top (Org, Contact Email, Exemption, Status, Keywords, Version, Contract#)
- Pipeline name with ph-core namespace
- Badges (CI/CD, nf-test, Nextflow version, nf-core template)
- Clear introduction describing pipeline purpose
- Workflow diagram or visualization
- Usage instructions with samplesheet format
- Parameters documentation
- Output description
- Credits and contributors
- Citations including nf-core framework
- DISCLAIMER reference

**meta.yml files must include:**
- Tools section with `args_id`
- Input/output types (map, file, directory, string, boolean, integer, float, list)
- Mandatory/Optional markings
- Java glob patterns for file types
- Proper descriptions
- Tools section listing every tool used in module
- Input/output entries matching corresponding channels in module
- Keywords sufficient for findability (research domain, data types, tool function)
- Input/output tuples split into separate entries
- Descriptive file content descriptions

## Review Process

### For Workflows
1. Check file structure and naming (lowercase, no punctuation)
2. Verify single-command execution
3. Review container references and versioning
4. Check documentation (README, CHANGELOG, citations)
5. Validate testing setup (nf-test, GitHub Actions)
6. Verify cloud compatibility configurations
7. Check branch structure
8. Review license and attribution
9. Verify SHARE IT Act metadata block
10. Validate workflow specificity (single pipeline per data/analysis type)
11. Check static assets and large test data stored separately
12. Verify no deprecated NextFlow features (entry points replaced with profiles)
13. Validate proper credits and acknowledgements
14. Check nf-core template version (v3.2.0 or later)
15. Verify custom branding (nf-core branding removed)
16. Validate Docker support with versioned software bundles
17. Check semantic versioning implementation
18. Verify standardized parameter usage
19. Validate minimum input requirements
20. Check presence of required directories (`assets/`, `bin/`, `conf/`, `docs/`, `libs/`)
21. Verify nf-core lint compliance
22. Validate nf-tools usage in GitHub Actions for version control
23. Check workflow runs on cloud environments (AWS/Azure/GCP)

### For Modules
1. Verify directory structure: `modules/*/main.nf` + `meta.yaml`
2. Check process naming (UPPERCASE with underscore)
3. Review container directive (specific version)
4. Validate input/output channel definitions
5. Check for stub block
6. Verify `versions.yml` emission
7. Review `$task.ext.args` usage
8. Check code alignment
9. Validate meta.yaml completeness
10. Review test coverage (nf-test)
11. Verify all mandatory input files in input channel definitions
12. Check all optional input files in input channel definitions
13. Validate non-mandatory arguments via `$task.ext.args` (not `val` channels)
14. Check `$args` variables correspond to tool order in command line
15. Verify `when` statement not changed in process definition
16. Check stub block exists for all modules
17. Verify no repeated `script:` blocks
18. Validate non-file mandatory commands as value channels
19. Check multi-command piping justified by runtime/storage advantage
20. Verify inputs/outputs compressed where applicable
21. Check exit code handling with `||` operator if tool doesn't exit with 0
22. Validate standard meta fields (`meta.id`, `meta.single_end`) or document non-standard fields
23. Verify versions.yml removes leading "v" from version numbers
24. Check module stored in `modules/` directory with proper naming
25. Validate container in AMDP Elastic Container Registry
26. Verify container passes Trivy security scans
27. Check only one process per module
28. Validate resource label appropriate (`process_single`, `process_low`, `process_medium`, `process_high`)
29. If multi-threading supported, verify `$task.cpus` usage
30. Check multi-tool containers available on BioContainers or justify new container

### For Subworkflows
1. Verify minimum 2 modules
2. Check directory naming (lowercase)
3. Review channel naming (`snake_case`, descriptive)
4. Validate input declarations in `take` block
5. Check output emissions (named, typed)
6. Review documentation and comments
7. Verify test coverage
8. Check input/output channel structure comments
9. Verify channel structure descriptions in meta.yml
10. Validate tags for dependent modules in tests
11. Check for `assertAll()` usage with success assertions
12. Verify named file extensions for outputs (unless dataset/directory)
13. Check optional input handling with `[]` workaround
14. Validate input channel names signify object type
15. Ensure output channel names based on major file type only
16. Verify no assumption of named params from parent workflow

### For README.md
1. Verify SHARE IT Act metadata block format and completeness
2. Check pipeline namespace (ph-core/*)
3. Validate badge links and accuracy
4. Ensure introduction is complete (no TODO comments)
5. Verify workflow diagram exists
6. Check usage section completeness (samplesheet format documented)
7. Validate parameters table
8. Check outputs documentation reference
9. Verify credits section is complete
10. Validate citations and nf-core attribution
11. Check DISCLAIMER reference

### For Test Data
1. Verify test data doesn't replicate existing test-data unnecessarily
2. Check test data is publicly available with appropriate licenses
3. Ensure test data is minimal size
4. Validate test data files described in branch's README (source, generation method, licenses)
5. Verify reuse of pre-existing files from nf-core/test-datasets where possible
6. Check test data files have entries in nf-core/test-datasets repo README
7. For modules requiring large data or local databases, verify stub feature usage
8. Validate files organized by discipline, organism, platform, or format
9. Check downstream/related test-data files named based on upstream file names

## Review Feedback Structure

When providing code review feedback, structure comments as:

1. **Critical Issues** (MUST fix - Tier 1 violations):
   - Specific file/line references
   - Explanation of why it's critical
   - Suggested fix

2. **Important Issues** (SHOULD fix - Tier 2 violations):
   - File/line references
   - Impact explanation
   - Improvement suggestions

3. **Recommendations** (NICE to have - Tier 3):
   - Best practice suggestions
   - Optimization opportunities

4. **Compliant Patterns**: Acknowledge what's done well

## Important Reminders

- Be specific with file names and line numbers
- Explain why something matters
- Provide code examples for fixes
- Check all applicable standards for the tier
- Acknowledge good patterns
- Focus on ph-core and AMDP compliance standards
- Ensure documentation completeness (no TODO placeholders in production code)