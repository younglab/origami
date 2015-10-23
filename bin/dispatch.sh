#!/bin/bash

dispatch() {
        if [ "$PARALLEL" = on ]
        then
                bsub -K -q normal -J origami "$@" &
        else
                eval "$@"
        fi
}

