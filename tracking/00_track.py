#!/usr/bin/env python3
"""
Project State Tracking System
Manages project_state.json for active task tracking
"""

import json
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent
STATE_FILE = PROJECT_ROOT / "project_state.json"


def load_state():
    """Load project state from JSON file"""
    if not STATE_FILE.exists():
        print(f"ERROR: State file not found: {STATE_FILE}")
        sys.exit(1)

    with open(STATE_FILE, 'r') as f:
        return json.load(f)


def save_state(state):
    """Save project state to JSON file"""
    state['last_updated'] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    with open(STATE_FILE, 'w') as f:
        json.dump(state, f, indent=2)

    print(f"✅ State saved: {state['last_updated']}")


def status(args=None):
    """Show current project status"""
    state = load_state()

    print("=" * 70)
    print(f"PROJECT: {state['project']}")
    print(f"LAST UPDATED: {state['last_updated']}")
    print(f"OVERALL PROGRESS: {state['overall_progress']}%")
    print("=" * 70)
    print()

    # Current phase
    current_phase_id = state['current_phase']
    current_phase = state['phases'][current_phase_id]
    print(f"🔄 CURRENT PHASE: {current_phase['name']}")
    print(f"   Status: {current_phase['status']}")
    print(f"   Progress: {current_phase['progress']}%")
    print(f"   {current_phase['description']}")
    print()

    # Show current phase tasks
    print("TASKS IN CURRENT PHASE:")
    for task in current_phase['tasks']:
        icon = "✅" if task['status'] == "complete" else "🔄" if task['status'] == "in_progress" else "⏳"
        print(f"  {icon} [{task['status']}] {task['name']}")
    print()

    # Blocking issues
    if state['blocking_issues']:
        print("🚫 BLOCKING ISSUES:")
        for issue in state['blocking_issues']:
            print(f"  - {issue}")
        print()

    # Reviewer requirements status
    print("REVIEWER REQUIREMENTS STATUS:")
    for reviewer, reqs in state['reviewer_requirements'].items():
        print(f"  {reviewer.upper()}:")
        for req, status in reqs.items():
            icon = "✅" if status == "complete" else "⏳"
            print(f"    {icon} {req}: {status}")
    print()

    # Phase summary
    print("ALL PHASES:")
    for phase_id, phase in state['phases'].items():
        icon = "✅" if phase['status'] == "complete" else "🔄" if phase['status'] == "in_progress" else "⏳"
        progress_bar = "█" * (phase['progress'] // 10) + "░" * (10 - phase['progress'] // 10)
        print(f"  {icon} {phase['name']}: [{progress_bar}] {phase['progress']}%")
    print()


def complete(args):
    """Mark a task as complete"""
    if not args:
        print("ERROR: Task ID required. Usage: track.py complete <task_id>")
        sys.exit(1)

    task_id = args[0]
    state = load_state()

    # Find and mark task complete
    task_found = False
    for phase_id, phase in state['phases'].items():
        for task in phase['tasks']:
            if task['id'] == task_id:
                task['status'] = "complete"
                task_found = True
                print(f"✅ Task marked complete: {task['name']}")

                # Update phase progress
                total_tasks = len(phase['tasks'])
                completed_tasks = sum(1 for t in phase['tasks'] if t['status'] == "complete")
                phase['progress'] = int((completed_tasks / total_tasks) * 100)

                # Check if phase is complete
                if completed_tasks == total_tasks:
                    phase['status'] = "complete"
                    print(f"🎉 Phase complete: {phase['name']}")

                    # Move to next phase
                    phase_ids = list(state['phases'].keys())
                    current_idx = phase_ids.index(phase_id)
                    if current_idx + 1 < len(phase_ids):
                        next_phase_id = phase_ids[current_idx + 1]
                        state['current_phase'] = next_phase_id
                        state['phases'][next_phase_id]['status'] = "in_progress"
                        print(f"📍 Moving to next phase: {state['phases'][next_phase_id]['name']}")

                # Update overall progress
                all_tasks = []
                for p in state['phases'].values():
                    all_tasks.extend(p['tasks'])
                completed = sum(1 for t in all_tasks if t['status'] == "complete")
                state['overall_progress'] = int((completed / len(all_tasks)) * 100)

                # Check reviewer requirements
                if 'reviewer_requirement' in task:
                    req_path = task['reviewer_requirement'].split('_')
                    reviewer = f"{req_path[0]}_{req_path[1]}_critical"
                    req_name = '_'.join(req_path[2:])
                    if reviewer in state['reviewer_requirements']:
                        state['reviewer_requirements'][reviewer][req_name] = "complete"
                        print(f"✅ Reviewer requirement satisfied: {task['reviewer_requirement']}")

                break
        if task_found:
            break

    if not task_found:
        print(f"ERROR: Task not found: {task_id}")
        print("Available task IDs:")
        for phase in state['phases'].values():
            for task in phase['tasks']:
                print(f"  - {task['id']}: {task['name']}")
        sys.exit(1)

    save_state(state)


def next_task(args=None):
    """Show next pending task"""
    state = load_state()
    current_phase = state['phases'][state['current_phase']]

    print("=" * 70)
    print("NEXT TASK TO WORK ON:")
    print("=" * 70)

    for task in current_phase['tasks']:
        if task['status'] == "pending":
            print(f"\n📌 Task ID: {task['id']}")
            print(f"📝 Task: {task['name']}")
            if 'reviewer_requirement' in task:
                print(f"📊 Reviewer Requirement: {task['reviewer_requirement']}")
            print(f"\n💡 To mark complete: ./tracking/00_track.py complete {task['id']}")
            return

    print("\n✅ All tasks in current phase complete!")
    print(f"Moving to next phase...")


def blockers(args=None):
    """Show blocking issues"""
    state = load_state()

    if not state['blocking_issues']:
        print("✅ No blocking issues")
        return

    print("🚫 BLOCKING ISSUES:")
    for issue in state['blocking_issues']:
        print(f"  - {issue}")


def report(args=None):
    """Generate progress report"""
    state = load_state()

    print("=" * 70)
    print(f"PROGRESS REPORT - {state['last_updated']}")
    print("=" * 70)
    print()

    print(f"Overall Progress: {state['overall_progress']}%")
    print()

    for phase_id, phase in state['phases'].items():
        total = len(phase['tasks'])
        completed = sum(1 for t in phase['tasks'] if t['status'] == "complete")
        in_progress = sum(1 for t in phase['tasks'] if t['status'] == "in_progress")
        pending = sum(1 for t in phase['tasks'] if t['status'] == "pending")

        print(f"{phase['name']}:")
        print(f"  Status: {phase['status']}")
        print(f"  Progress: {phase['progress']}%")
        print(f"  Tasks: {completed}/{total} complete, {in_progress} in progress, {pending} pending")
        print()

    print("Reviewer Requirements:")
    for reviewer, reqs in state['reviewer_requirements'].items():
        complete_count = sum(1 for status in reqs.values() if status == "complete")
        total_count = len(reqs)
        print(f"  {reviewer}: {complete_count}/{total_count} complete")
    print()


def main():
    if len(sys.argv) < 2:
        print("Usage: track.py <command> [args]")
        print("\nCommands:")
        print("  status       - Show current project status")
        print("  complete <task_id> - Mark task as complete")
        print("  next         - Show next pending task")
        print("  blockers     - Show blocking issues")
        print("  report       - Generate progress report")
        sys.exit(1)

    command = sys.argv[1]
    args = sys.argv[2:] if len(sys.argv) > 2 else None

    commands = {
        'status': status,
        'complete': complete,
        'next': next_task,
        'blockers': blockers,
        'report': report
    }

    if command not in commands:
        print(f"ERROR: Unknown command: {command}")
        sys.exit(1)

    commands[command](args)


if __name__ == "__main__":
    main()
