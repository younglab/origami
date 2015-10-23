#!/bin/bash

BINDIR=~/dsday/origami/bin ### Need to generalize this
OUTPUTDIR=output
VERBOSE=off
SKIP=on
PARALLEL=off
BZPOSTFIX="[.]bz2$"

source $BINDIR/dispatch.sh

verbose() {
	if [ "$VERBOSE" = on ]
	then
		NOWTIME=$(date)
		echo "[$NOWTIME] $1"
	fi
}

TEMP=`getopt -o o::hvap -l output::,noskip -n 'origami' -- "$@"`
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
		--noskip)
			SKIP=off
			;;
		-p)
			PARALLEL=on
			;;
	esac
	shift
done

LEFTREADS="$1"
RIGHTREADS="$2"

echo "Launching origami..."

verbose "Analyzing $LEFTREADS and $RIGHTREADS"

verbose "Creating output directory"
mkdir $OUTPUTDIR
verbose "Creating temporary file directory"
mkdir $OUTPUTDIR/tmp

### handle zip status
if [[ $LEFTREADS =~ $BZPOSTFIX ]]
then
	dispatch "bzcat $LEFTREADS > $OUTPUTDIR/tmp/left_unzip.fq"
	LEFTREADS=$OUTPUTDIR/tmp/left_unzip.fq
fi

if [[ $RIGHTREADS =~ $BZPOSTFIX ]]
then
        dispatch "bzcat $RIGHTREADS > $OUTPUTDIR/tmp/right_unzip.fq"
        RIGHTREADS=$OUTPUTDIR/tmp/right_unzip.fq
fi

wait

verbose "Removing adapter sequences on $LEFTREADS and $RIGHTREADS"
[ "$SKIP" = off -o ! -e "$OUTPUTDIR/tmp/left_kept.fq" ] && $BINDIR/adapter_trim.sh $OUTPUTDIR/tmp $PARALLEL $LEFTREADS $RIGHTREADS


verbose "Aligning reads"
[ "$SKIP" = off -o ! -e "$OUTPUTDIR/tmp/left_kept.bam" ] && $BINDIR/bowtie_align.sh $OUTPUTDIR $PARALLEL
echo "Done"
