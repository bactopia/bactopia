#!/usr/bin/env bash
set -euo pipefail

# Wrapper script to locate and run bactopia-status.
# Search order:
#   1. bactopia-status on PATH
#   2. conda env named "bactopia-py"
#   3. conda env matching "bactopia-py*" (first match)
#   4. Fail with helpful error

find_conda_env() {
    local envs
    envs=$(conda env list --json 2>/dev/null \
        | python3 -c "import sys,json; [print(p) for p in json.load(sys.stdin)['envs']]" 2>/dev/null) || return 1

    # Prefer exact "bactopia-py" env
    local exact
    exact=$(echo "$envs" | grep '/bactopia-py$' | head -1)
    if [[ -n "$exact" ]]; then
        echo "$exact"
        return 0
    fi

    # Fall back to any bactopia-py* env
    local fuzzy
    fuzzy=$(echo "$envs" | grep '/bactopia-py' | head -1)
    if [[ -n "$fuzzy" ]]; then
        echo "$fuzzy"
        return 0
    fi

    return 1
}

# 1. Check PATH
if command -v bactopia-status &>/dev/null; then
    exec bactopia-status "$@"
fi

# 2. Search conda envs
if env_path=$(find_conda_env); then
    exec conda run --prefix "$env_path" bactopia-status "$@"
fi

# 3. Give up
echo "ERROR: bactopia-status not found." >&2
echo "Install bactopia-py or activate a conda env containing it." >&2
exit 1
