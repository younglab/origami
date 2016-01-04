#!/bin/bash

dispatch() {
        if [ "$PARALLEL" = on ]
        then
                bsub -K -q normal -J origami -o logs/cluster_log.txt "$@" &
        else
                eval "$@"
        fi
}

