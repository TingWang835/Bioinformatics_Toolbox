### Note: this is prepared for possible future AI agent implication (with openclaw in mind)
## Directory Overview
- **./Snakefile**: Master control; handles project switching.
- **./toolbox/**: Modular `.smk` logic files (Code/Rules).
- **./refs/[ref_name]/**: Shared references and indexes (Read-only for most rules).
- **./reads/[project_name]/**: Individual project workspace.
    - **config.yaml**: Project parameters.
    - **[results]/**: Sub-folders for BAM, VCF, and logs.

## Maintenance & AI Protocol
- **Logging**: All rules must use `log:` with `2>&1` redirection.
- **Pathing**: Use paths relative to the root directory to ensure portability.
- **Reference Integrity**: Do not allow rules to modify files in `./refs/` unless specifically designed for indexing.

## Execution Logic
1. **Toolbox Isolation**: Rules in `toolbox/` should use relative paths based on the root `Snakefile`.
2. **Data Specificity**: Always verify the `config.yaml` path inside the `Snakefile` before running a new project.
3. **Outputs**: All results must be directed to the specific project subdirectory within `reads/` to avoid cross-contamination.
