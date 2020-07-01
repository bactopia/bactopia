#!/bin/bash
set -e
set -u

mkdir amrdb/ amrdb-temp/
amrfinder_update -d amrdb-temp/
find amrdb-temp/ -name "database_format_version.txt" | xargs -I {} dirname {} | xargs -I {} sh -c 'mv {}/* amrdb/'
tar -czvf amrdb.tar.gz amrdb/
