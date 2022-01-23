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
    mpirun -np $(($i+1)) ./SubGorder $GRAPH $SAMPLE_RATIO $i $RATIO $AVERAGE
    for j in `seq 1 $AVERAGE`
    do
      ./../ligra/utils/SNAPtoAdj -s $i"-"$j"-sample.txt" $i"-"$j".adj"
      ./../ligra/apps/PageRank $i"-"$j".adj" | grep -o "[0-9.]*" >> $i."exectime"
    done
    rm *sample.txt
    rm *.adj
    if [ $i -eq 1 ]
    then
      echo "Original"
      ./../ligra/utils/SNAPtoAdj "original.txt" "original.adj"
      ./../ligra/apps/PageRank "original.adj" | grep -o "[0-9.]*" >> "original.exectime"
      rm original.txt
      rm original.adj
      #echo "Sequential"
      #./../ligra/utils/SNAPtoAdj "seq.txt" "seq.adj"
      #./../ligra/apps/PageRank "seq.adj" | grep -o "[0-9.]*" > "seq.exectime"
      rm seq.txt
      rm seq.adj
      #echo "Random"
      #./../ligra/utils/SNAPtoAdj "random.txt" "random.adj"
      #./../ligra/apps/PageRank "random.adj" | grep -o
      #rm random.txt
      #rm random.adj
    fi
    #mpirun -np $(($i+1)) ./Gorder-025 $GRAPH $SAMPLE_RATIO $i $AVERAGE
    #mpirun -np $(($i+1)) ./Gorder-050 $GRAPH $SAMPLE_RATIO $i $AVERAGE
    #mpirun -np $(($i+1)) ./Gorder-075 $GRAPH $SAMPLE_RATIO $i $AVERAGE
  done
  #python3 evaluation.py $PARALLEL
  rm *.txt
  rm *.ans
  #rm *.exectime
fi
