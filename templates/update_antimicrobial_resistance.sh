#!/bin/bash
set -e
set -u

AMRFINDER_DB=`which amrfinder | sed 's=bin/amrfinder=share/amrfinderplus/data='`
if [ ! -d "$AMRFINDER_DB" ]; then
    # DB does not exist, create it
    amrfinder -u
else
    if [ "!{params.update_amr}" == "true" ]; then
        # User asked for update
        amrfinder -u
    fi
fi
