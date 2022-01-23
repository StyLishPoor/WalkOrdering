#!/bin/bash
CORE=`lscpu | grep socket | grep -o "[0-9]*"`
HYPER=`lscpu | grep core | grep -o "[0-9]*"`
MAX_THREAD=`expr $(($CORE*$HYPER))`
SELECT_GRAPH=$1
SAMPLE_RATIO=$2
PARALLEL=$3
MINIMUM_RATIO=$4
AVERAGE=$5
rm *.txt
rm *.ans
rm *.exectime
rm *.rwtime
if [ "$(($3+1))" -gt "$MAX_THREAD" ]
then
  echo "Exceed Max Thread $MAX_THREAD"
else
  echo "My Proposal"
  mpirun -np $(($3+1)) ./SubGorder-eval-final $SELECT_GRAPH $SAMPLE_RATIO $3 $MINIMUM_RATIO $AVERAGE
  for j in `seq 1 $AVERAGE`
  do
    ./../ligra/utils/SNAPtoAdj $3"-"$j"-sample.txt" $3"-"$j".adj"
    #./../ligra/apps/PageRank $3"-"$j".adj" | grep -o "[0-9.]*" >> $TRUE_GRAPH"-"$i."exectime"
    ./../ligra/apps/PageRank $3"-"$j".adj" | grep -o "[0-9.]*"
    ./../ligra/apps/PageRankDelta $3"-"$j".adj" | grep -o "[0-9.]*"
  done
  rm *.adj
  rm *.txt
  rm *.ans
  rm *.exectime
  rm *.rwtime
  echo "Before"
  mpirun -np $(($3+1)) ./SubGorder-eval-final $SELECT_GRAPH $SAMPLE_RATIO $3 1 $AVERAGE
  for j in `seq 1 $AVERAGE`
  do
    ./../ligra/utils/SNAPtoAdj $3"-"$j"-sample.txt" $3"-"$j".adj"
    #./../ligra/apps/PageRank $3"-"$j".adj" | grep -o "[0-9.]*" >> $TRUE_GRAPH"-"$i."exectime"
    ./../ligra/apps/PageRank $3"-"$j".adj" | grep -o "[0-9.]*"
    ./../ligra/apps/PageRankDelta $3"-"$j".adj" | grep -o "[0-9.]*"
  done
  rm *.adj
  rm *.txt
  rm *.txt
  rm *.ans
  rm *.exectime
  rm *.rwtime
fi
