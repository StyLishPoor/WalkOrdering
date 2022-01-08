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
    echo "----[$i]----"
    mpirun -np $(($i+1)) ./Gorder $GRAPH $SAMPLE_RATIO $i $AVERAGE
    less $i"-sample.txt" | wc -l
    ./../ligra/utils/SNAPtoAdj $i"-sample.txt" $i".adj"
    ./../ligra/apps/PageRank $i".adj"
    #./../ligra/apps/PageRankDelta $i".adj"
    #./../ligra/apps/BC $i".adj"
    #./../ligra/apps/CF $i".adj"
    #./../ligra/apps/KCore $i".adj"
    rm $i"-sample.txt"
    rm $i".adj"
    if [ $i -eq 1 ]
    then
      echo "Original"
      ./../ligra/utils/SNAPtoAdj "original.txt" "original.adj"
      ./../ligra/apps/PageRank "original.adj"
      #./../ligra/apps/PageRankDelta "original.adj"
      #./../ligra/apps/BC "original.adj"
      #./../ligra/apps/CF "original.adj"
      #./../ligra/apps/KCore "original.adj"
      rm original.txt
      rm original.adj
      echo "Sequential"
      ./../ligra/utils/SNAPtoAdj "seq.txt" "seq.adj"
      ./../ligra/apps/PageRank "seq.adj"
      #./../ligra/apps/PageRankDelta "seq.adj"
      #./../ligra/apps/BC "seq.adj"
      #./../ligra/apps/CF "seq.adj"
      #./../ligra/apps/KCore "seq.adj"
      rm seq.txt
      rm seq.adj
      echo "Random"
      ./../ligra/utils/SNAPtoAdj "random.txt" "random.adj"
      ./../ligra/apps/PageRank "random.adj"
      #./../ligra/apps/PageRankDelta "random.adj"
      #./../ligra/apps/BC "random.adj"
      #./../ligra/apps/CF "random.adj"
      #./../ligra/apps/KCore "random.adj"
      rm random.txt
      rm random.adj
    fi
    #mpirun -np $(($i+1)) ./Gorder-025 $GRAPH $SAMPLE_RATIO $i $AVERAGE
    #mpirun -np $(($i+1)) ./Gorder-050 $GRAPH $SAMPLE_RATIO $i $AVERAGE
    #mpirun -np $(($i+1)) ./Gorder-075 $GRAPH $SAMPLE_RATIO $i $AVERAGE
  done
  #python3 evaluation.py $PARALLEL
  rm *.ans
fi
