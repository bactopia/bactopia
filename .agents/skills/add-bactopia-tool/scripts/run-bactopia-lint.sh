#!/usr/bin/env bash
set -euo pipefail

# Wrapper script to run bactopia-lint scoped to a single tool across all three tiers.
# Usage: run-bactopia-lint.sh <tool_name> --bactopia-path <path> [--pretty]
#
# Runs lint on:
#   --module <tool>                        (module tier)
#   --subworkflow <tool>                   (subworkflow tier)
#   --workflow bactopia-tools/<tool>       (workflow tier)
#
# Search order for bactopia-lint:
#   1. bactopia-lint on PATH
#   2. conda env named "bactopia-dev"
#   3. conda env named "bactopia-py"
#   4. conda env matching "bactopia-*" (first match)
#   5. Fail with helpful error

find_conda_env() {
    local envs
    envs=$(conda env list --json 2>/dev/null \
        | python3 -c "import sys,json; [print(p) for p in json.load(sys.stdin)['envs']]" 2>/dev/null) || return 1

    local dev_exact
    dev_exact=$(echo "$envs" | grep '/bactopia-dev$' | head -1)
    if [[ -n "$dev_exact" ]]; then
        echo "$dev_exact"
        return 0
    fi

    local exact
    exact=$(echo "$envs" | grep '/bactopia-py$' | head -1)
    if [[ -n "$exact" ]]; then
        echo "$exact"
        return 0
    fi

    local fuzzy
    fuzzy=$(echo "$envs" | grep '/bactopia-' | head -1)
    if [[ -n "$fuzzy" ]]; then
        echo "$fuzzy"
        return 0
    fi

    return 1
}

resolve_lint() {
    if command -v bactopia-lint &>/dev/null; then
        echo "bactopia-lint"
        return 0
    fi

    if env_path=$(find_conda_env); then
        echo "conda run --prefix $env_path bactopia-lint"
        return 0
    fi

    echo "ERROR: bactopia-lint not found." >&2
    echo "Install bactopia-py or activate a conda env containing it." >&2
    exit 1
}

# Parse arguments
TOOL=""
BACTOPIA_PATH=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bactopia-path)
            BACTOPIA_PATH="$2"
            shift 2
            ;;
        --pretty|--json|--quiet)
            EXTRA_ARGS+=("$1")
            shift
            ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            TOOL="$1"
            shift
            ;;
    esac
done

if [[ -z "$TOOL" ]]; then
    echo "Usage: run-bactopia-lint.sh <tool_name> --bactopia-path <path> [--pretty]" >&2
    exit 1
fi

if [[ -z "$BACTOPIA_PATH" ]]; then
    echo "ERROR: --bactopia-path is required" >&2
    exit 1
fi

LINT_CMD=$(resolve_lint)

# Run lint scoped to all three tiers for this tool
# shellcheck disable=SC2086
$LINT_CMD \
    --module "$TOOL" \
    --subworkflow "$TOOL" \
    --workflow "bactopia-tools/$TOOL" \
    --bactopia-path "$BACTOPIA_PATH" \
    "${EXTRA_ARGS[@]}"
