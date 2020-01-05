for folder in *
do
 bsub -q 18 /lab/solexa_bartel/usage/find_all_large_txt_in_a_folder.sh $folder 1000000 ./usage/$folder-large-files-20180123
done
