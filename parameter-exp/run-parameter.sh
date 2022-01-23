#!/bin/bash
CORE=`lscpu | grep socket | grep -o "[0-9]*"`
HYPER=`lscpu | grep core | grep -o "[0-9]*"`
MAX_THREAD=`expr $(($CORE*$HYPER))`
SAMPLE_RATIO=$1
AVERAGE=$2
FLAG=0
rm *.txt
rm *.ans
rm *.exectime
rm *.rwtime
if [ "$(($3+1))" -gt "$MAX_THREAD" ]
then
  echo "Exceed Max Thread $MAX_THREAD"
else
  for GRAPH in ~/data/evaluation/*.txt
  do
    TMP_GRAPH=${GRAPH: 0:-4}
    TRUE_GRAPH=${TMP_GRAPH##*/}
    echo $TRUE_GRAPH
    for i in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
    #for i in `seq 1 $PARALLEL`
    do 
      echo "----[$i]----"
      mpirun -np 6 ./SubGorder-eval-parameter $GRAPH $SAMPLE_RATIO $i $AVERAGE
      for j in `seq 1 $AVERAGE`
      do
        ~/ligra/utils/SNAPtoAdj $i"-"$j"-sample.txt" $i"-"$j".adj"
        ~/ligra/apps/PageRank $i"-"$j".adj" | grep -o "[0-9.]*" >> $TRUE_GRAPH"-"$i."exectime"
      done
      rm *sample.txt
      rm *.adj
      if [ $FLAG -eq 0 ]
      then
        FLAG=1
        echo "Original"
        ~/ligra/utils/SNAPtoAdj "original.txt" "original.adj"
        ~/ligra/apps/PageRank "original.adj" | grep -o "[0-9.]*" >> $TRUE_GRAPH"-original.exectime"
        echo "Random"
        ~/ligra/utils/SNAPtoAdj "random.txt" "random.adj"
        ~/ligra/apps/PageRank "random.adj" | grep -o "[0-9.]*" >> $TRUE_GRAPH"-random.exectime"
        rm original.txt
        rm original.adj
        rm random.txt
        rm random.adj
      fi
    done
    FLAG=0
  done
  python3 evaluation-parameter.py
  #rm *.txt
  rm *.ans
  rm *.exectime
  rm *.rwtime
fi
