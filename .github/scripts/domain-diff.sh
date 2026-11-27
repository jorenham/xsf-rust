#!/usr/bin/env bash
# Runs domain_diff.py in CI and writes its table to the job summary.
#
# Compares the merge commit with the base tip for the functions that depend on the changed files
# of src/. Exit code 1 is a regression. Exit code 2 is a comparison that could not be made (e.g. a
# new function, which the base lacks, or a missing reference), which only warns. Exit code 3 is a
# function that got slower, which fails unless the PR has the label `allow-slowdown`; re-run the
# job after adding or removing it (a label event starts no run). The last line of the table tells
# these from a failure of the script itself (e.g. a build error).
#
# Needs PR (the number of the PR) and GH_TOKEN and GH_REPO for `gh`.
set -euo pipefail

table=$RUNNER_TEMP/table.md
code=0
uv run --project tools tools/accuracy/domain_diff.py \
  --base HEAD^1 --affected --max-inputs 2097152 > "$table" || code=$?
cat "$table"
{
  echo "## Accuracy over whole domains"
  cat "$table"
} >> "$GITHUB_STEP_SUMMARY"
case "$code:$(tail -n 1 "$table")" in
  0:*) ;;
  1:*FAIL)
    echo "::error::A function got less accurate over a domain (see the job summary and the artifact domain-diff)"
    exit 1
    ;;
  2:*INCONCLUSIVE) echo "::warning::The comparison is incomplete (see the job summary)" ;;
  3:*SLOW)
    # read live: the labels of the event are those from when the run started
    query='any(.labels[]; .name == "allow-slowdown")'
    if ! allowed=$(gh pr view "$PR" --json labels --jq "$query"); then
      echo "::error::The labels of the PR could not be read"
      exit 1
    elif [ "$allowed" != true ]; then
      echo "::error::A function got slower (see the job summary); add the label allow-slowdown and re-run the job if that is intended"
      exit 1
    fi
    echo "::warning::A function got slower (see the job summary); allowed by the label allow-slowdown"
    ;;
  *)
    echo "::error::domain_diff.py failed with exit code $code"
    exit 1
    ;;
esac
