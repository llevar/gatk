#!/usr/bin/env bash

. utils.sh

RESULTS_CSV=test_case_md_results.csv

for exec_cores in 1; do
  for total_cores in 64; do
    num_exec=$((total_cores * exec_cores))
    time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded --shardedOutput true" $num_exec $exec_cores 4g 4g >> $RESULTS_CSV
    time_gatk "SortReadFileSpark -I hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/q4_spark_eval/out/markdups" $num_exec $exec_cores 4g 4g
  done
done