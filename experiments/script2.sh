#!/bin/bash
# Experimentos con INF-EFPS

TIMEOUT=900
NRUN=5
MAXNODES=15000
declare -a SOLVERS=("efpss")
declare -a CUTMAX=("5" "10" "25")
declare -a CUTFREQ=("100" "1000" "10000")
CASES="instances.csv"
BIN="../pds-lim"
INPUT="../inputs/"
OUTPUT="../outputs/exp2/"

mkdir -p $OUTPUT 

tail -n +2 $CASES | while read -r line 
do
    IFS=',' read -r index name vertices degree <<< $line
    if [ $vertices -gt $MAXNODES ]
    then
        continue
    fi
    # solve
    for solver in "${SOLVERS[@]}"
    do
        for omega in $(seq 0 $degree)
        do
            date
            time $BIN -s $solver -w $omega -f $INPUT$name -n $NRUN -t $TIMEOUT -o $OUTPUT --in-prop
            time $BIN -s $solver -w $omega -f $INPUT$name -n $NRUN -t $TIMEOUT -o $OUTPUT --out-prop
            time $BIN -s $solver -w $omega -f $INPUT$name -n $NRUN -t $TIMEOUT -o $OUTPUT --in-prop --out-prop
            time $BIN -s $solver -w $omega -f $INPUT$name -n $NRUN -t $TIMEOUT -o $OUTPUT --out-prop --init-efps
            time $BIN -s $solver -w $omega -f $INPUT$name -n $NRUN -t $TIMEOUT -o $OUTPUT --in-prop --out-prop --init-efps
            for cutmax in "${CUTMAX[@]}"
            do
                for cutfreq in "${CUTFREQ[@]}"
                do
                    time $BIN -s $solver -w $omega -f $INPUT$name -n $NRUN -t $TIMEOUT -o $OUTPUT --out-prop --init-efps --cuts --cut-max $cutmax --cut-freq $cutfreq
                done
            done
        done
    done
done
