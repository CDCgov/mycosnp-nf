# GitHub Copilot Custom Instructions

This directory contains custom instructions for GitHub Copilot to provide context-specific guidance for code reviews and development in this repository.

## Structure

- `*.instructions.md` - Custom instruction files that guide Copilot's behavior

## Current Instructions

### ph-core-review.instructions.md

**Purpose:** Provides ph-core and AMDP compliance standards for code reviews

**Applies to:** All files in the repository (`applyTo: "**"`)

**Excluded agents:** Coding agent (only used by code review)

**Key features:**

- Three-tier compliance checking (AMDP Prerequisites, ph-core Requirements, ph-core Recommendations)
- NextFlow DSL2 standards enforcement
- Module, subworkflow, and workflow compliance validation
- Container security and versioning requirements (Trivy scans)
- Documentation standards (README with SHARE IT Act metadata, CHANGELOG, meta.yml)
- Testing requirements (nf-test framework with stub tests)
- Detailed review checklists (23 workflow items, 30 module items, 16 subworkflow items)
- Test data management standards
- Code alignment standards ('Harshil Alignment™️')
- Resource labeling requirements
- nf-core template compliance (v3.2.0+)

## Usage Contexts

These instructions are used by GitHub Copilot in:

1. **VS Code Editor**: When chatting with Copilot or requesting code reviews locally
2. **GitHub.com**: When Copilot performs automated pull request reviews
3. **Code Review Agent**: Applies standards during PR code review process

Note: Currently, there's no way to enable instructions only for GitHub.com while disabling them in VS Code - the instructions apply in both contexts.

## References

- [GitHub Copilot Custom Instructions Documentation](https://docs.github.com/copilot/how-tos/configure-custom-instructions/add-repository-instructions)
- [Agent-Specific Instructions Announcement](https://github.blog/changelog/2025-11-12-copilot-code-review-and-coding-agent-now-support-agent-specific-instructions/)