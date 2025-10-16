# Session Context - October 16, 2025

## Current Status
- Working on QC extraction from nf-core outputs
- Just cloned collaborator repo to data/collaborator_repos/
- Created comparison document between our approach and theirs
- Ready to run QC analysis scripts

## Key Files Created/Updated Today
1. `project_status.md` - Updated with current progress
2. `analysis_strategy.md` - Comprehensive plan for the project
3. `analysis_comparison.md` - Comparison with collaborator methods
4. `scripts/reorganization_plan.md` - Plan to restructure scripts
5. `results/01_qc_analysis/` - Directory structure created

## Next Immediate Steps
1. Run QC extraction: `./scripts/03_qc_analysis/run_qc_analysis.sh`
2. Review extracted metrics
3. Compare with collaborator QC standards

## Important Paths
- Container: `/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif`
- Outputs: `output/PRJNA*/SRR*/`
- Collaborator repo: `data/collaborator_repos/`

## Container Commands
```bash
# Always start with:
module load singularity/3.9.3
cd /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq
singularity shell --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif
cd /proj
```

## Active Todo Items
- Extract QC metrics from all 44 samples (HIGH PRIORITY)
- Create QC summary tables and figures
- Compare our approach with collaborators' methods (IN PROGRESS)
- Import 34 GWAS candidate genes list
- Reorganize scripts directory

## Questions/Decisions Pending
- Which reference genome version to use for final analysis
- How to handle platform batch effects (GAIIx vs HiSeq2000)
- Where to get the 34 GWAS candidate genes list