#!/usr/bin/env python3
"""
Visual Progress Display for albopictus-diapause-rnaseq
Creates a clear visual representation of project progress
"""

import json
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).parent.parent
STATE_FILE = PROJECT_ROOT / "tracking" / "project_state.json"

def load_state():
    with open(STATE_FILE, 'r') as f:
        return json.load(f)

def create_visual():
    state = load_state()

    output = []

    # Header
    output.append("="*80)
    output.append("ALBOPICTUS DIAPAUSE RNA-SEQ PROGRESS DASHBOARD")
    output.append("="*80)
    output.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    output.append(f"Overall Progress: {state['overall_progress']}%")
    output.append("="*80)

    # THE THREE-PART ANALYSIS STRATEGY
    output.append("\n📋 THE ANALYSIS STRATEGY (3 PARTS):")
    output.append("-" * 80)
    output.append("PART 1: REPLICATE COLLABORATORS (featureCounts)")
    output.append("        → Validate pipeline on new AalbF3 genome")
    output.append("        → Use collaborator's exact DESeq2 scripts")
    output.append("")
    output.append("PART 1.5: VALIDATE SALMON VS FEATURECOUNTS")
    output.append("        → Must prove r > 0.95 correlation")
    output.append("        → If fails: stop, do not proceed to Part 2")
    output.append("")
    output.append("PART 2: REVIEWER REQUIREMENTS (Salmon - only after validation)")
    output.append("        → Genome-wide statistics")
    output.append("        → Fisher's exact test for GWAS enrichment")
    output.append("")
    output.append("PART 3: GWAS CANDIDATES (last)")
    output.append("        → Extract candidate expression")
    output.append("        → Only after Parts 1 & 2 complete")

    # Current Work
    output.append("\n🎯 CURRENT TASK:")
    output.append("-" * 80)
    current_phase = state['phases'][state['current_phase']]
    found_task = False
    for task in current_phase['tasks']:
        if task['status'] in ['pending', 'failed']:
            output.append(f"➡️  {task['name']}")
            output.append(f"    ID: {task['id']}")
            if 'script' in task:
                output.append(f"    📄 Script: {task['script']}")
            if 'output_directory' in task:
                output.append(f"    📁 Output: {task['output_directory']}")
            if 'note' in task:
                output.append(f"    ⚠️  Note: {task['note']}")

            # Show exact command to run
            if 'script' in task:
                script_path = task['script']
                if script_path.endswith('.sh'):
                    output.append(f"    🔧 Run: sbatch {script_path}")
                elif script_path.endswith('.py'):
                    output.append(f"    🔧 Run: singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c \"cd /proj && python {script_path}\"")
                elif script_path.endswith('.R'):
                    output.append(f"    🔧 Run: singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c \"cd /proj && Rscript {script_path}\"")

            output.append(f"    ✅ Mark complete: python tracking/00_track.py complete {task['id']}")
            found_task = True
            break

    if not found_task:
        output.append("✅ All tasks in current phase complete!")
        output.append("    Check status: python tracking/00_track.py status")

    # Phase Progress Bars
    output.append("\n📊 PHASE PROGRESS:")
    output.append("-" * 80)

    for phase_id, phase in state['phases'].items():
        # Count tasks
        total = len(phase['tasks'])
        complete = sum(1 for t in phase['tasks'] if t['status'] == 'complete')
        failed = sum(1 for t in phase['tasks'] if t['status'] == 'failed')

        # Create progress bar
        bar_width = 30
        filled = int((complete / total) * bar_width) if total > 0 else 0
        bar = "█" * filled + "░" * (bar_width - filled)

        # Status indicator
        if phase['status'] == 'complete':
            indicator = "✅"
        elif phase['status'] == 'failed':
            indicator = "❌"
        elif phase['status'] == 'in_progress':
            indicator = "🔄"
        else:
            indicator = "⏳"

        # Format phase name to fixed width
        phase_name = phase['name'][:50].ljust(50)

        # Print phase line
        output.append(f"{indicator} {phase_name} [{bar}] {complete}/{total}")

        # Show current task if in progress or failed
        if phase['status'] in ['in_progress', 'failed']:
            for task in phase['tasks']:
                if task['status'] in ['pending', 'in_progress', 'failed']:
                    status_icon = "⚠️" if task['status'] == 'failed' else "🔄"
                    task_line = f"    {status_icon} {task['name']}"
                    if 'script' in task:
                        task_line += f" ({task['script']})"
                    output.append(task_line)
                    break

    # Reviewer Requirements
    output.append("\n🎓 REVIEWER REQUIREMENTS:")
    output.append("-" * 80)
    req = state['reviewer_requirements']

    output.append("Reviewer 2 (CRITICAL):")
    output.append(f"  • Genome-wide stats: {req['reviewer2_critical']['genome_wide_stats']}")
    output.append(f"  • Fisher's exact test: {req['reviewer2_critical']['fishers_exact_test']}")

    output.append("\nReviewer 3 (CRITICAL):")
    output.append(f"  • Sample table: {req['reviewer3_critical']['sample_table']}")
    output.append(f"  • DESeq2 design explanation: {req['reviewer3_critical']['deseq2_design_explanation']}")
    output.append(f"  • Stage-specific expression: {req['reviewer3_critical']['stage_specific_expression']}")
    output.append(f"  • Mapping statistics: {req['reviewer3_critical']['mapping_statistics']}")

    # Key Files
    output.append("\n📁 KEY FILES:")
    output.append("-" * 80)
    output.append("• Pipeline output: output/star_salmon/featurecounts/*.tsv (22,177 genes)")
    output.append("• Container: albopictus-diapause-rnaseq.sif (1.5GB, has all packages)")
    output.append("• Tracking: tracking/00_track.py + tracking/project_state.json")
    output.append("• Project info: tracking/project_summary.md")

    # Footer
    output.append("\n" + "="*80)
    output.append("NEXT STEP: Run 'tracking/00_track.py next' to see current task")
    output.append("="*80)

    # Write to file
    output_file = PROJECT_ROOT / 'tracking' / 'visual_progress.txt'
    with open(output_file, 'w') as f:
        f.write('\n'.join(output))

    # Also print to console
    print('\n'.join(output))
    print(f"\n✅ Visual progress saved to: {output_file}")

if __name__ == "__main__":
    create_visual()