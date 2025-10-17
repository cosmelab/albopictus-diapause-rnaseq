# Git Hook Setup - AI Attribution Blocker

## What Was Done

Created a global git commit-msg hook that automatically rejects any commit messages containing AI attribution (Claude, Anthropic, etc.).

## Problem Solved

Prevents accidental commits that include:
- "Claude"
- "Anthropic"
- "Co-Authored-By: Claude"
- "Generated with Claude"
- AI assistant references
- 🤖 emoji

## Implementation Steps

### 1. Created Global Template Directory
```bash
mkdir -p ~/.git-templates/hooks
```

### 2. Created commit-msg Hook
```bash
cat > ~/.git-templates/hooks/commit-msg << 'EOF'
#!/bin/bash
# Reject commits with AI attribution

commit_msg_file="$1"
commit_msg=$(cat "$commit_msg_file")

# Check for forbidden terms
if echo "$commit_msg" | grep -iE "(Claude|Anthropic|Co-Authored-By: Claude|Generated with.*Claude|🤖|AI assistant)" > /dev/null; then
    echo "ERROR: Commit message contains AI attribution!"
    echo "This violates repository rules."
    echo "Please remove any references to Claude, Anthropic, or AI assistance."
    exit 1
fi

exit 0
EOF
```

### 3. Made Hook Executable
```bash
chmod +x ~/.git-templates/hooks/commit-msg
```

### 4. Configured Git Globally
```bash
git config --global init.templateDir ~/.git-templates
```

### 5. Applied to Current Repository
```bash
git init
```

## How It Works

- **New repositories**: Automatically get the hook when created or cloned
- **Existing repositories**: Run `git init` in the repo directory to apply the template
- **Current repository (ucr-ento-social)**: Already protected

## Testing the Hook

Try to commit with AI attribution:
```bash
git commit -m "Test commit

Co-Authored-By: Claude <noreply@anthropic.com>"
```

Expected result:
```
ERROR: Commit message contains AI attribution!
This violates repository rules.
Please remove any references to Claude, Anthropic, or AI assistance.
```

## Applying to Other Existing Repositories

Navigate to any existing repository and run:
```bash
cd /path/to/your/repo
git init
```

This will copy the hook from the template directory without affecting your existing git history.

## Location of Files

- **Global template**: `~/.git-templates/hooks/commit-msg`
- **Repository hook**: `.git/hooks/commit-msg` (in each repo)
- **Git config**: Run `git config --global --get init.templateDir` to verify

## Modifying the Hook

To update the hook for all future repositories:
1. Edit `~/.git-templates/hooks/commit-msg`
2. Run `git init` in existing repositories to update them

## Reference

This hook enforces rules from `assistant_rules.md`:
- Rule 2: Never include AI assistant names, "Anthropic", or any AI attribution in commits
- Rule 3: Never use commit messages that mention AI assistance
