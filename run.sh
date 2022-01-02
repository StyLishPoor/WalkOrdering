#!/bin/bash
CORE=`lscpu | grep socket | grep -o "[0-9]*"`
HYPER=`lscpu | grep core | grep -o "[0-9]*"`
MAX_THREAD=`expr $(($CORE*$HYPER))`
GRAPH=$1
SAMPLE_RATIO=$2
PARALLEL=$3
AVERAGE=$4
if [ "$(($3+1))" -gt "$MAX_THREAD" ]
then
  echo "Exceed Max Thread $MAX_THREAD"
else
  for i in `seq 1 $PARALLEL`
  do 
    mpirun -np $(($i+1)) ./Gorder $GRAPH $SAMPLE_RATIO $i $AVERAGE
    mpirun -np $(($i+1)) ./Gorder-025 $GRAPH $SAMPLE_RATIO $i $AVERAGE
    mpirun -np $(($i+1)) ./Gorder-050 $GRAPH $SAMPLE_RATIO $i $AVERAGE
    mpirun -np $(($i+1)) ./Gorder-075 $GRAPH $SAMPLE_RATIO $i $AVERAGE
  done
  python3 evaluation.py $PARALLEL
  rm *.ans
fi
