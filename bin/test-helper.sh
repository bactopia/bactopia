#! /bin/bash
set -euo pipefail
set -x

# first pass
BACTOPIA_TESTS=/home/rpetit3/repos/bactopia/bactopia-tests/ nf-test test main.nf.test

# second pass
BACTOPIA_TESTS=/home/rpetit3/repos/bactopia/bactopia-tests/ nf-test test main.nf.test

# cleanup
rm -rf .nf-test/ .nf-test.log
