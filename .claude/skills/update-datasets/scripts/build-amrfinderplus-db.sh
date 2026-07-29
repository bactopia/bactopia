#!/usr/bin/env bash
set -euo pipefail

# Build the AMRFinder+ database tarball inside the amrfinderplus module's pinned
# container. This script owns the build recipe formerly carried by the (now
# removed) modules/amrfinderplus/update/ Nextflow module.
#
# It reads the pinned container image + tool version from the run module's
# config (modules/amrfinderplus/run/module.config) so a tool version bump flows
# through automatically, then runs `amrfinder_update` in that container and
# packages the result as amrfinderplus.tar.gz.

usage() {
    cat >&2 <<'EOF'
Usage: build-amrfinderplus-db.sh --bactopia-path <path> [--runtime docker|singularity] [--outdir <dir>]

  --bactopia-path   Path to the Bactopia repository (required).
  --runtime         Container runtime: docker (default) or singularity.
  --outdir          Output directory for the tarball
                    (default: <bactopia-path>/logs/dataset-builds/<timestamp>).
EOF
    exit 1
}

BACTOPIA_PATH=""
RUNTIME="docker"
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bactopia-path) BACTOPIA_PATH="${2:-}"; shift 2 ;;
        --runtime)       RUNTIME="${2:-}"; shift 2 ;;
        --outdir)        OUTDIR="${2:-}"; shift 2 ;;
        -h|--help)       usage ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage ;;
    esac
done

[[ -n "$BACTOPIA_PATH" ]] || { echo "ERROR: --bactopia-path is required" >&2; usage; }
[[ -d "$BACTOPIA_PATH" ]] || { echo "ERROR: bactopia path not found: $BACTOPIA_PATH" >&2; exit 1; }

case "$RUNTIME" in
    docker|singularity) ;;
    *) echo "ERROR: --runtime must be 'docker' or 'singularity'" >&2; exit 1 ;;
esac

MODULE_CONFIG="$BACTOPIA_PATH/modules/amrfinderplus/run/module.config"
[[ -f "$MODULE_CONFIG" ]] || { echo "ERROR: module.config not found: $MODULE_CONFIG" >&2; exit 1; }

# Pinned images live in the run module's config (ext.docker / ext.image).
DOCKER_IMAGE=$(grep -oP 'ext\.docker\s*=\s*"\K[^"]+' "$MODULE_CONFIG" | head -1)
SING_IMAGE=$(grep -oP 'ext\.image\s*=\s*"\K[^"]+' "$MODULE_CONFIG" | head -1)
[[ -n "$DOCKER_IMAGE" ]] || { echo "ERROR: could not parse ext.docker from $MODULE_CONFIG" >&2; exit 1; }
[[ -n "$SING_IMAGE" ]]   || { echo "ERROR: could not parse ext.image from $MODULE_CONFIG" >&2; exit 1; }

# ext.docker is a bare Docker Hub-style repo (e.g. biocontainers/...); Nextflow
# prepends params.registry (default quay.io) at run time, so we do the same. If
# the image's first path segment already looks like a host (contains '.' or ':'),
# it is left untouched.
REGISTRY=$(grep -oP 'registry\s*=\s*"\K[^"]+' "$BACTOPIA_PATH/conf/params.config" | head -1)
REGISTRY="${REGISTRY:-quay.io}"
first_segment="${DOCKER_IMAGE%%/*}"
if [[ "$DOCKER_IMAGE" == "$first_segment" || ( "$first_segment" != *.* && "$first_segment" != *:* ) ]]; then
    FULL_IMAGE="$REGISTRY/$DOCKER_IMAGE"
else
    FULL_IMAGE="$DOCKER_IMAGE"
fi

if [[ -z "$OUTDIR" ]]; then
    TS=$(date +%Y%m%d_%H%M%S)
    OUTDIR="$BACTOPIA_PATH/logs/dataset-builds/$TS"
fi
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)

# Recipe executed inside the container; /work is bound to $OUTDIR. Command
# substitutions ($()) are evaluated by the container's shell at runtime.
read -r -d '' BODY <<'EOF' || true
set -euo pipefail
cd /work
rm -rf amrfinderplus-temp amrfinderplus
mkdir -p amrfinderplus-temp
amrfinder_update -d amrfinderplus-temp
mv "amrfinderplus-temp/$(readlink amrfinderplus-temp/latest)" amrfinderplus/
tar czvf amrfinderplus.tar.gz amrfinderplus/
amrfinder --version > TOOL_VERSION
echo $(amrfinder --database amrfinderplus --database_version 2> /dev/null) | rev | cut -f1 -d' ' | rev > DB_VERSION
rm -rf amrfinderplus-temp amrfinderplus
EOF

echo ">> Building AMRFinder+ database" >&2
echo ">> runtime=$RUNTIME image=$FULL_IMAGE outdir=$OUTDIR" >&2

if [[ "$RUNTIME" == "docker" ]]; then
    command -v docker >/dev/null 2>&1 || { echo "ERROR: docker not found" >&2; exit 1; }
    docker run --rm --user "$(id -u):$(id -g)" \
        -v "$OUTDIR":/work -w /work "$FULL_IMAGE" \
        bash -c "$BODY"
else
    RUNNER=""
    if command -v apptainer >/dev/null 2>&1; then
        RUNNER="apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        RUNNER="singularity"
    else
        echo "ERROR: singularity/apptainer not found" >&2
        exit 1
    fi
    # Reuse the same registry-qualified image via docker:// so both runtimes match.
    "$RUNNER" exec --bind "$OUTDIR":/work "docker://$FULL_IMAGE" bash -c "$BODY"
fi

TARBALL="$OUTDIR/amrfinderplus.tar.gz"
[[ -s "$TARBALL" ]] || { echo "ERROR: build produced no tarball at $TARBALL" >&2; exit 1; }

SHA=$(sha256sum "$TARBALL" | cut -d' ' -f1)
TOOL_VERSION=$(cat "$OUTDIR/TOOL_VERSION" 2>/dev/null || echo "unknown")
DB_VERSION=$(cat "$OUTDIR/DB_VERSION" 2>/dev/null || echo "unknown")

cat <<EOF

=== AMRFinder+ database build summary ===
tarball=$TARBALL
sha256=$SHA
tool_version=$TOOL_VERSION
db_version=$DB_VERSION
image=$FULL_IMAGE
runtime=$RUNTIME
EOF
