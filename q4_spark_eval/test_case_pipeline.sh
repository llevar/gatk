#!/usr/bin/env bash

. utils.sh

RESULTS_CSV=test_case_pipeline_results.csv

run_id=0
for exec_cores in 1; do
  for total_cores in 64; do
    num_exec=$((total_cores * exec_cores))
    #time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded --shardedOutput true" $num_exec $exec_cores 4g 4g >> $RESULTS_CSV
    time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded --shardedOutput true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --knownSites hdfs:///user/$USER/q4_spark_eval/dbsnp_138.b37.vcf --joinStrategy SHUFFLE" $num_exec $exec_cores 4g 4g >> $RESULTS_CSV
    #time_gatk "ReadsPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/reads-sharded --shardedOutput true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --knownSites hdfs:///user/$USER/q4_spark_eval/dbsnp_138.vcf" $num_exec $exec_cores 4g 4g >> $RESULTS_CSV
    run_id=$((run_id+1))
  done
done