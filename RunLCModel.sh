RunLCModelOn=$1
tmp_dir=$2

CurrentComputer=$(hostname)
if ! [[ $CurrentComputer == $RunLCModelOn ]]; then
	echo -e "\nThere seems to be a problem to ssh to computer $RunLCModelOn."
	echo -e"Did you forget to take your meds? Or did you forget to provide a key, so that ssh is possible without password authentication?"
fi


echo -e "\n\n7. LCmodel processing started!\n\n"
chmod 755 ${tmp_dir}/lcm_*
KillLCMFile="${tmp_dir}/KillLCMProcesses.sh"
echo "echo 'Kill lcm-processes.'" > $KillLCMFile 
chmod 775 $KillLCMFile
pidloopy=0
for process_file in ${tmp_dir}/lcm_process_core_*; do
	echo $process_file
	let pidloopy=pidloopy+1;
	$process_file &
	pid_list[pidloopy]=$!
	echo -e "kill ${pid_list[pidloopy]}" >> $KillLCMFile 
done
wait  #for all cores to finish their jobs

unset pid_list;
rm $KillLCMFile
