#! /bin/bash

find $1 -type f | sort | xargs -I {} md5sum {} | awk '{print "    - path: "$2"\n      md5sum: "$1}'
