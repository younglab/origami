#!/bin/bash

BINDIR=~/dsday/origami/bin/ ### Need to generalize this
OUTPUTDIR=output

TEMP=`getopt -o o::h -l output:: -n 'origami' -- "$@"`
eval set -- "$TEMP"

while [ $# -ge 1 ]; do
	case "$1" in
		--)
			shift
			break
			;;
		-o|--output)
			OUTPUTDIR=$2
			shift
			;;
		-h)
			echo "Help menu"
			exit 0
			;;
	esac
	shift
done

echo "Launching origami..."

mkdir $OUTPUTDIR
mkdir $OUTPUTDIR/tmp


$BINDIR/adapter_trim.sh $OUTPUTDIR/tmp $1 $2
