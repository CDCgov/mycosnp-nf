#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
script=$SCRIPT_DIR/mycosnp_create_sample_sheet.pl
echo "sample,r1,r2,r3,r4"
for VAR in "$@"
do
    $script -i $VAR -f
    for DIR in `ls $VAR`; do
        if [[ -d "$VAR/$DIR" ]]; then
            $script -i $VAR/$DIR -f
        fi
    done
done | sort