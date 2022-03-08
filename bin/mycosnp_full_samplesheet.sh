#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
script=$SCRIPT_DIR/mycosnp_create_sample_sheet.pl
echo "sample,r1,r2,r3,r4"
$script -i $1 -f
for DIR in `ls $1`; do
    if [[ -d "$1/$DIR" ]]; then
        $script -i $1/$DIR -f
    fi
done | sort 