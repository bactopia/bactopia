#!/usr/bin/env bash
set -euo pipefail

# Wrapper script to locate and run bactopia-catalog.
# Search order:
#   1. bactopia-catalog on PATH
#   2. conda env named "bactopia-dev" (preferred -- built from environment.yml at the bactopia repo root)
#   3. conda env named "bactopia-py"
#   4. conda env matching "bactopia-*" (first match)
#   5. Fail with helpful error

find_conda_env() {
    local envs
    envs=$(conda env list --json 2>/dev/null \
        | python3 -c "import sys,json; [print(p) for p in json.load(sys.stdin)['envs']]" 2>/dev/null) || return 1

    # Prefer exact "bactopia-dev" env (built from environment.yml at the bactopia repo root)
    local dev_exact
    dev_exact=$(echo "$envs" | grep '/bactopia-dev$' | head -1)
    if [[ -n "$dev_exact" ]]; then
        echo "$dev_exact"
        return 0
    fi

    # Fall back to exact "bactopia-py" env
    local exact
    exact=$(echo "$envs" | grep '/bactopia-py$' | head -1)
    if [[ -n "$exact" ]]; then
        echo "$exact"
        return 0
    fi

    # Fall back to any bactopia-* env (catches bactopia-py-dev, bactopia-dev-v4, etc.)
    local fuzzy
    fuzzy=$(echo "$envs" | grep '/bactopia-' | head -1)
    if [[ -n "$fuzzy" ]]; then
        echo "$fuzzy"
        return 0
    fi

    return 1
}

# Default --llms-output so llms.txt is always regenerated alongside catalog.json
# (bactopia-catalog only renders llms.txt when the flag is given).
if [[ " $* " != *" --llms-output "* ]]; then
    llms_target="llms.txt"
    prev=""
    for arg in "$@"; do
        case "$prev" in --bactopia-path) llms_target="$arg/llms.txt"; break ;; esac
        case "$arg" in --bactopia-path=*) llms_target="${arg#*=}/llms.txt"; break ;; esac
        prev="$arg"
    done
    set -- "$@" --llms-output "$llms_target"
fi

# 1. Check PATH
if command -v bactopia-catalog &>/dev/null; then
    exec bactopia-catalog "$@"
fi

# 2. Search conda envs
if env_path=$(find_conda_env); then
    exec conda run --prefix "$env_path" bactopia-catalog "$@"
fi

# 3. Give up
echo "ERROR: bactopia-catalog not found." >&2
echo "Install bactopia-py or activate a conda env containing it." >&2
exit 1
