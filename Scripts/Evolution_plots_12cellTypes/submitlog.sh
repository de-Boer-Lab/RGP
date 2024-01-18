#!/bin/bash -l
#insert use commands here
#  $1 is the task file


export id=`awk "NR==$PBS_ARRAY_INDEX" $1` # input the $PBS_ARRAY_INDEXth line of this file into $id
#uncomment and modify the next two lines if your task file has multiple fields.
splitID=($id) #split id on whitespace into array that can be accessed like ${splitID[0]}
export track_name=${splitID[0]} # if there was whitesapce delim fields in the task file, make the new task ID the first field (0)
export cell_type=${splitID[1]}
export enformer_track=${splitID[2]}
export track_index=${splitID[3]}
export bed_file=${splitID[5]}
echo  $PBS_ARRAY_INDEX : $id #print the current task to the primary output log
#scratch only
export logDir="/scratch/st-cdeboer-1/iluthra/randomDNA/Analysis_outputs_06152023/H1_DNase_rep1_07052023/Evolution_plots_12cellTypes/log_dir/" #fill this in with where you would like the output logs to go - this directory must actually exist.
export logPre=$logDir/$track_index #a prefix to use for the job output log ($logPre.olog) and the file denoting that the job completed ($logPre.done)
#redirect stdout and stderr to file
exec 1>>$logPre.olog # redirect standard out to an output log file based off of this task ID
exec 2>&1 #redirect stderr to stdout
set -e #make it so that the script exits completely if one command fails

module load gcc python miniconda3 cuda cudnn
conda activate /arc/project/st-cdeboer-1/rafi_tf
conda activate ishika_enformer

prediction_path="/scratch/st-cdeboer-1/iluthra/randomDNA/all_iPSC_related_H1rep1/"
output_path_dir="/scratch/st-cdeboer-1/iluthra/randomDNA/Analysis_outputs_06152023/H1_DNase_rep1_07052023/"
scripts_path="/project/st-cdeboer-1/iluthra/enformer/test_enformer/de-former/Paper_Review_Analysis/Final_clean_scripts_figures/Evolution_plots_12cellTypes/"


if [ ! -e $logPre.done ] #makes sure the stuff after "then" is only run if the job was not completed last time.
then
        echo Job $PBS_JOBID:$PBS_ARRAY_INDEX started on $PBS_O_HOST: `date` #print to the output file the current job id and task, and host
        #PLACE COMMANDS HERE
        #echo python "$scripts_path""encode_enformer_bed_files_allCells.py" --path_to_predictions "$prediction_path" --cell_type $cell_type --track_index $track_index --track_name $track_name --encode_file_paths $bed_file --output_path "$output_path_dir""/Evolution_plots_12cellTypes/"

        python "$scripts_path""encode_enformer_bed_files_allCells.py" --path_to_predictions "$prediction_path" --cell_type $cell_type --track_index $track_index --track_name $track_name --encode_file_paths $bed_file --output_path "$output_path_dir""/Evolution_plots_12cellTypes/"

        conda deactivate
        module load gcc/5.5.0
        module load bedtools2/2.30.0
        cd "$output_path_dir""/Evolution_plots_12cellTypes/"

        bedtools intersect -wao -a genomic_regions_coordinates_n1000_114kb_128bins_"$cell_type"_"$track_name".bed -b "$cell_type"_"$track_name"_encode.bed > "$cell_type"_"$track_name"_enformer_encode_overlap.bed

        conda deactivate
        module load gcc python miniconda3 cuda cudnn
        conda activate ishika_enformer

        python "$scripts_path""Evolution_plots_all_celltypes_final.py" --path_to_predictions "$prediction_path" --cell_type $cell_type --track_index $track_index --track_name $track_name --bed_file "$cell_type"_"$track_name"_enformer_encode_overlap.bed --output_path "$output_path_dir""/Evolution_plots_12cellTypes/" > 80_percent_out_"$cell_type"_"$track_name".txt

        touch $logPre.done #always the last command - create the file that indicates this task completed successfully
        echo Job finished: `date`
else
        echo Job already done. To redo, rm $logPre.done
fi
qstat -r $PBS_JOBID | tail -n 1
