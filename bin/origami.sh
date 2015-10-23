#!/bin/bash

BINDIR=~/dsday/origami/bin ### Need to generalize this
OUTPUTDIR=output
VERBOSE=off
SKIP=on

verbose() {
	if [ "$VERBOSE" = on ]
	then
		NOWTIME=$(date)
		echo "[$NOWTIME] $1"
	fi
}

TEMP=`getopt -o o::hva -l output:: -n 'origami' -- "$@"`
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
		-v)
			VERBOSE=on
			;;
		-a)
			SKIP=off
			;;
	esac
	shift
done

echo "Launching origami..."

verbose "Creating output directory"
mkdir $OUTPUTDIR
verbose "Creating temporary file directory"
mkdir $OUTPUTDIR/tmp


verbose "Removing adapter sequences"
[ "$SKIP" = off -o ! -e "$OUTPUTDIR/tmp/left_kept.fq" ] && $BINDIR/adapter_trim.sh $OUTPUTDIR/tmp $1 $2


verbose "Aligning reads"
[ "$SKIP" = off -o ! -e "$OUTPUTDIR/tmp/left_kept.bam" ] && $BINDIR/bowtie_align.sh $OUTPUTDIR
echo "Done"
