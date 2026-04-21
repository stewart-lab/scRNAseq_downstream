#!/bin/bash

# run_cellchat_parallel.sh
# Usage: ./run_cellchat_parallel.sh <filelist> <data_dir> <group_by> <source1> <source2> <start_task_id> <end_task_id>
# Example: ./run_cellchat_parallel.sh filelist2.txt /w5home/bmoore/Pierre_sc_zebrafish/ celltype 1 2 1 5

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <filelist> <data_dir> <group_by> <source1> <source2> <start_task_id> <end_task_id>"
    exit 1
fi

FILELIST=$1
DATA_DIR=$2
GROUP_BY=$3
SOURCE1=$4
SOURCE2=$5
START=$6
END=$7

for i in $(seq $START $END); do
    echo "Submitting task $i for filelist $FILELIST"
    nohup Rscript cellchat.R --filelist "$FILELIST" --data_dir "$DATA_DIR" --group_by "$GROUP_BY" --source1 "$SOURCE1" --source2 "$SOURCE2" --task_id $i > "cellchat_task_${i}.log" 2>&1 &
done

echo "Submitted tasks $START to $END. Check logs cellchat_task_*.log for progress."
