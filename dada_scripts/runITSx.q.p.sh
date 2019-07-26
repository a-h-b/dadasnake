#! /bin/bash

usage() {
    echo "Usage: $0 [-o <string>] [-t <string>] [-f <file suffix>] [-i <ITSregion>] [-n <number of ITSxregions>] [-e <e-value cut-off>] [-d] [-r] [-g] [-h] input " 1>&2
    echo "       input should be one of the following: " 1>&2 
    echo "       i) an absolute path to a fasta file " 1>&2
    echo "       ii)a file containing the absolute path(s) to the fasta file(s)"  1>&2
    echo "       iii) an absolute path to the directory with the fasta file(s)" 1>&2
    echo "       -o <string> the string is the name of a directory that will be made to save  all ITSx outputs; defaults to <ITSx_output>"  1>&2
    echo "       -t <string> directory where all output goes or it is  written to the place where the script was called" 1>&2
    echo "       -f <file suffix> if a directory is specified as the input, only files with this suffix will be used; defaults to fasta" 1>&2
    echo "       -i <ITSregion> ITS1 or ITS2; defaults to ITS2" 1>&2
    echo "       -n <number of ITSregions> 1-4; defaults to 1" 1>&2
    echo "       -e <e-value cutoff> e-value; defaults to 1e-5" 1>&2
    echo "       -d if set, ITSx will give more detailed output" 1>&2
    echo "       -r if set, all ITSx output will be removed and only the pruned fungal sequences will be kept; not recommended" 1>&2
    echo "       -g if set, (and -r isn't) the ITSx output will also contain a graphical representation of all identified ITS2 sequences" 1>&2
    echo "       -h displays this message" 1>&2 
    echo " " 1>&2
    exit 1 
}

while getopts o:t:rgdf:i:n:e:h flag
do
    case $flag in
        o)
            OUTPUT_STR=$OPTARG;;
        t)
            OUTDIR_STR=$OPTARG;;
        r)  
	    REMOVE_TMP=T;;
	d)
	    DETAILS=T;;
	g)
            GRAPH=T;;
	f)
	    SUFF=$OPTARG;;
	i)
	    ITS_REG=$OPTARG;;
        n)
            ITS_NUM=$OPTARG;;
        e)
            EVALUE=$OPTARG;;
        h)
            usage
	    exit ;;
        :) 
	    echo "Missing option argument for -$OPTARG" >&2 
	    usage
	    exit 1;;
        *)  
	    echo "Unimplemented option: -$OPTARG" >&2 
            usage
	    exit 1;;
        ?)
	    usage
            exit
	     ;;
    esac
done

shift $((OPTIND-1))

if [ -z "$OUTPUT_STR" ]
then
    OUTPUT_STR="ITSx_out"
fi

if [ -z "$OUTDIR_STR" ]
then
    OUTDIR_STR="."
fi

if [ -z "$ITS_REG" ]
then
    ITS_REG="ITS2"
fi

if [ -z "$ITS_NUM" ]
then
    ITS_NUM=1
fi

if [ -z "$EVALUE" ]
then
    EVALUE="1e-5"
fi

if [ -z "$DETAILS" ]
then
    DETAILS=F
fi

if [ -z "$REMOVE_TMP" ]
then
    REMOVE_TMP=F
fi

if [ -z "$GRAPH" ]
then
    GRAPH=F
fi

if [ -z "$1" ];
then
    echo "missing input"
    usage
    exit 1
else
    INFILE=$1
fi


#if the file cannot be found
if [[ !  -e "$1" ]]; then
   echo "Inputfile "$1" was not found."
   echo "Provide full path of input file or input directory"
   echo "or try qsub -cwd <Skript> <File>"
   exit 1
fi

# check if there are multiple files given
if [[  "$#" -gt 1 ]]; then
   echo "Number of input files is $#"
   echo "This program only takes one input file but offers options to process multiple sequence files"
   usage
   exit 1
fi


if [[ -d $INFILE ]]
then
    if [ -z ${SUFF+x} ]
    then
        SUFF="fasta"
    elif [ -z "$SUFF" ]
    then
        echo "Suffix cannot be an empty string, if target files are specified by directory"
        exit 1
    fi
    CURRDIR=$(pwd)
    cd $INFILE
    INPUTS=(*.$SUFF)
    cd $CURRDIR
    for i in ${!INPUTS[@]}
    do
	INPUTS[i]="${INFILE}/${INPUTS[i]}"
    done
elif [[ -f $INFILE ]]
then
    if ! [ -z ${SUFF+x} ]
    then
	echo "Warning: Suffix was given but will be ignored"
	unset SUFF
    fi

# check if the input file contains a list of full paths or is a fasta file
# if input is a fasta file, get the full path and write that to a default file which is used as input
    if [[ $(cat $INFILE | grep -E ^\/) || $(head -n 1 $INFILE | grep -E ^\>)  ]]; then
       if [[ $(head -n 1 $INFILE | grep -E ^\>) ]]; then
          abspath=$(cd "$(dirname "$1")" &>/dev/null  &&  printf "%s/%s" "$PWD" "${1##*/}") 
          INPUTS=($abspath) 
          #INFILE=ITS_input 
       else
	  readarray -t INPUTS < $INFILE
	fi

       SUFFA=()
       for i in ${INPUTS[@]}
       do
           currsuff="${i##*.}"
           SUFFA+=$currsuff
        done
        tmpsuff=($(echo "${SUFFA[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
        if [ ${#tmpsuff[@]} -eq 1 ]
        then
	     SUFF=${SUFFA[0]}
             unset SUFFA
        fi
    else
	echo "no proper input file." 
        exit 1
    fi   
else
    echo "invalid input."
    exit 1
fi

if [ ${#INPUTS[@]} -eq 0 ] || [ "${INPUTS[0]}" = "*.$SUFF" ]
then
    echo "no files selected"
    exit
else

module load hmmer/3.1b2-1
CURRDIR=$(pwd)
mkdir $OUTDIR_STR/$OUTPUT_STR
cd $OUTDIR_STR/$OUTPUT_STR
for i in "${!INPUTS[@]}"
do  
    FASTAFILE=${INPUTS[$i]}
    if [ -f $FASTAFILE ]
    then
	echo "$FASTAFILE"
	tmpout=${FASTAFILE##*/}
	if ! [ -z "$SUFFA" ]
	then
	    SUFF=${SUFFA[$i]}
	fi
        OUTNAME=${tmpout%.$SUFF}
	CMD="/data/project/metaamp/TOOLS/ITSx_1.0.11/ITSx -i $FASTAFILE --cpu ${NSLOTS:-1} --detailed_results $DETAILS --save_regions $ITS_REG --graphical $GRAPH -o $OUTNAME -N $ITS_NUM -E $EVALUE"
	echo $CMD
	eval $CMD 
	CMD="grep '|F|$ITS_REG' -A 1 --no-group-separator $OUTNAME.$ITS_REG.fasta | sed 's/|.*//' > ../$OUTNAME.$ITS_REG.fungi.fasta"
	echo $CMD
	eval $CMD
    else
	echo "Warning: Missing input file $FASTAFILE."
    fi
done

cd ..
if [ $REMOVE_TMP = "T" ]
then
    rmdir $OUTPUT_STR
echo "DONE"

fi

fi

