#!/bin/bash
CORE=`lscpu | grep socket | grep -o "[0-9]*"`
HYPER=`lscpu | grep core | grep -o "[0-9]*"`
MAX_THREAD=`expr $(($CORE*$HYPER))`
GRAPH=$1
SAMPLE_RATIO=$2
PARALLEL=$3
RATIO=$4
AVERAGE=$5
if [ "$(($3+1))" -gt "$MAX_THREAD" ]
then
  echo "Exceed Max Thread $MAX_THREAD"
else
  for i in `seq 1 $PARALLEL`
  do 
    echo "----[$i]----"
    mpirun -np $(($i+1)) ./SubGorder-eval $GRAPH $SAMPLE_RATIO $i $RATIO $AVERAGE
    for j in `seq 1 $AVERAGE`
    do
      ./../ligra/utils/SNAPtoAdj $i"-"$j"-sample.txt" $i"-"$j".adj"
      ./../ligra/apps/PageRank $i"-"$j".adj" | grep -o "[0-9.]*" >> $i."exectime"
    done
    rm *sample.txt
    rm *.adj
    if [ $i -eq 1 ]
    then
      echo "Original"
      ./../ligra/utils/SNAPtoAdj "original.txt" "original.adj"
      ./../ligra/apps/PageRank "original.adj" | grep -o "[0-9.]*" >> "original.exectime"
      echo "Random"
      ./../ligra/utils/SNAPtoAdj "random.txt" "random.adj"
      ./../ligra/apps/PageRank "random.adj" | grep -o "[0-9.]*" >> "random.exectime"
      rm original.txt
      rm original.adj
      rm random.txt
      rm random.adj
    fi
  done
  python3 evaluation.py $PARALLEL
  rm *.txt
  rm *.ans
  rm *.exectime
fi
