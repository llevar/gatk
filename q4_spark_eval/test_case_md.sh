#!/usr/bin/env bash

time_gatk() {
  GATK_ARGS=$1
  NUM_EXECUTORS=$2
  EXECUTOR_CORES=$3
  EXECUTOR_MEMORY=$4
  DRIVER_MEMORY=$5
  COMMAND=$(echo $GATK_ARGS | awk '{print $1}')
  LOG=${COMMAND}_$(date +%Y%m%d_%H%M%S).log
  ./gatk-launch $GATK_ARGS \
    -- \
    --sparkRunner SPARK --sparkMaster yarn-client --sparkSubmitCommand spark2-submit \
    --num-executors $NUM_EXECUTORS --executor-cores $EXECUTOR_CORES --executor-memory $EXECUTOR_MEMORY \
    --driver-memory $DRIVER_MEMORY \
  > $LOG 2>&1
  RC=$?
  DURATION_MINS=$(grep 'Elapsed time' $LOG | grep -Eow "[0-9]+\.[0-9][0-9]")
  echo "$GATK_ARGS,$NUM_EXECUTORS,$EXECUTOR_CORES,$EXECUTOR_MEMORY,$RC,$DURATION_MINS"
}

for num_exec in 64
do
  time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/markdups" $num_exec 1 4g 4g >> test_case_md_results.csv
done