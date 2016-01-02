#!/bin/bash

BINDIR=~/dsday/origami/bin ### Need to generalize this
OUTPUTDIR=output
VERBOSE=off
SKIP=on
PARALLEL=off
SPLITNUM=4000000
BZPOSTFIX="[.]bz2$"
BOWTIEIDX=notafile

source $BINDIR/dispatch.sh

verbose() {
	if [ "$VERBOSE" = on ]
	then
		NOWTIME=$(date)
		echo "[$NOWTIME] $1"
	fi
}

TEMP=`getopt -o o::hvap -l output::,noskip,splitnum::,bowtieidx: -n 'origami' -- "$@"`
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
		--splitnum)
		  SPLITNUM=$(expr "$2" \* 4)
		  shift
		  ;;
		--bowtieidx)
		  BOWTIEIDX=$2
		  shift
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
verbose "Creating logs directory"
mkdir $OUTPUTDIR/logs

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
[ "$SKIP" = off -o ! -e "$OUTPUTDIR/mapped_reads.bam" ] && $BINDIR/adapter_trim.sh $OUTPUTDIR $PARALLEL $SPLITNUM $LEFTREADS $RIGHTREADS

rm -f $OUTPUTDIR/tmp/left_unzip.fq  $OUTPUTDIR/tmp/right_unzip.fq

verbose "Aligning reads"
[ "$SKIP" = off -o ! -e "$OUTPUTDIR/mapped_reads.bam" ] && $BINDIR/bowtie_align.sh $OUTPUTDIR $BOWTIEIDX $PARALLEL $SPLITNUM

wait #finish all remaining processes

echo "Calling peaks"
$BINDIR/peak-calling.sh $OUTPUTDIR

#echo "Finding links"

#bedtools pairtobed -bedpe -type both -a $OUTPUTDIR/mapped_reads.bam -b $OUTPUTDIR/peaks_peaks.narrowPeak > $OUTPUTDIR/raw-links.out

echo "Done"
