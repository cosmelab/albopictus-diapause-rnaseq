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

## GitHub Authentication
```bash
# Load saved GitHub token
export GITHUB_TOKEN=$(cat ~/.github_token)
# Token saved securely at ~/.github_token (created Oct 16, 2025)
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

## TODO for QC Analysis (Next Session):
1. Keep PRJNA IDs but add clear labels in figures:
   - PRJNA268379 = "Huang et al. 2015 - Adult Females"
   - PRJNA158021 = "Poelchau et al. 2013a - Embryos"
   - PRJNA187045 = "Poelchau et al. 2013b - Pharate Larvae"

2. For QC, analyze per project first:
   - QC metrics table for PRJNA268379 (16 samples)
   - QC metrics table for PRJNA158021 (12 samples)
   - QC metrics table for PRJNA187045 (17 samples)

3. Create project-specific QC figures:
   - PCA for each project separately (to check for outliers within experiment)
   - Box plots of QC metrics grouped by experimental conditions
   - For PRJNA268379: group by photoperiod and blood meal
   - For PRJNA158021: group by timepoint and diapause status
   - For PRJNA187045: group by timepoint and diapause status

4. Then one overall figure showing all projects for comparison