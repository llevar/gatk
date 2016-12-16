#!/usr/bin/env bash

time_gatk() {
  GATK_ARGS=$1
  NUM_EXECUTORS=$2
  EXECUTOR_CORES=$3
  EXECUTOR_MEMORY=$4
  COMMAND=$(echo $GATK_ARGS | awk '{print $1}')
  LOG=${COMMAND}_$(date +%Y%m%d_%H%M%S).log
  ./gatk-launch $GATK_ARGS \
    -- \
    --sparkRunner SPARK --sparkMaster yarn-client --sparkSubmitCommand spark2-submit \
    --num-executors $NUM_EXECUTORS --executor-cores $EXECUTOR_CORES --executor-memory $EXECUTOR_MEMORY \
  > $LOG 2>&1
  RC=$?
  DURATION_MINS=$(grep 'Elapsed time' $LOG | grep -Eow "[0-9]+\.[0-9][0-9]")
  echo "$GATK_ARGS,$NUM_EXECUTORS,$EXECUTOR_CORES,$EXECUTOR_MEMORY,$RC,$DURATION_MINS"
}

echo 'Command,Number of executors,Executor cores,Executor memory,Exit code,Time (mins)' > test_case_4_results.csv

for num_exec in 4 8 16 32 64
do
  time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -L Broad.human.exome.b37.interval_list" $num_exec 1 1g >> test_case_4_results.csv
done

for num_exec in 2 4 8 16 32
do
  time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -L Broad.human.exome.b37.interval_list" $num_exec 2 1g >> test_case_4_results.csv
done

for num_exec in 1 2 4 8 16
do
  time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -L Broad.human.exome.b37.interval_list" $num_exec 4 1g >> test_case_4_results.csv
done