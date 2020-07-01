#!/bin/bash
set -e
set -u

mkdir amrdb/
amrfinder_update -d amrdb/
tar -czvf amrdb.tar.gz amrdb/
