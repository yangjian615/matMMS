#! /bin/bash

#
# Process a list of burst or survey files.
#
#   Calling Sequence
#     bash mms_fsm_l3_process_list.sh list
#     bash mms_fsm_l3_process_list.sh list log_file
#
# Inputs
#   LIST:       in, required
#               A file name with a list of burst or survey files to process. The
#                 file should have three columns: spacecraft, data rate mode, and
#                 start time of the data file.
#   LOG_FILE:   in, optional
#               The name of a file that will act as the batch processing log file.
#
# Exit Status
#   Upon exit, the exit status can take on a range of values from 0 to 255. Zero (0)
#     indicates no error, the values 1-99 indicate a warning, and 100-255 indicate an
#     error. Below are the specific error codes:
#         0  -  Everything OK
#       100  -  Bad input parameters
#       101  -  File not found
#       255  -  Unknown error not caught by MATLAB script
#

# Check inputs
if [ $# -ne 1 ] && [ $# -ne 2 ]; then
	echo 'Incorrect number of arguments.' 1>&2
	exit 100
fi

# Parse the argument list
list=${1}

# Log file
if [ $# -eq 2 ]; then
	logfile=${2}

	# Header for the log file
	cat > $logfile <<EOF
=========================
Log file for:       mms_fsm_l3_sdc
LOG_PATH_ROOT:      $LOG_PATH_ROOT
Log Filename:       $logfile
Runtime:            `date -u +%Y%m%d_%H%M%S`
Pwd:                `pwd`
Host Information:   `uname -a`
Arguments:          $0 $@
=========================
Calling MATLAB
EOF
else
	logfile=/dev/null
fi

# Loop over each line in LIST
while read -r sc mode tstart; do
	# Write the command to the log file
	echo "mms_fsm_l3_sdc( '$sc', '$mode', '$tstart' );" >> $logfile
	
	# Run the command
	#   - Must take input from null or it will consume the input file
	#   - Direct output to log file
	matlab -r -nosplash -nodesktop "\"try status = mms_fsm_l3_sdc( '$sc', '$mode', '$tstart' ); catch; exit(255); end; exit(status)\"" < /dev/null >> $logfile
	exitcode=$?
	
	# Indicate an error has occurred.
	if [ ${exitcode} -gt 100 ]; then
		echo "ERROR: $exitcode -- $sc $mode $tstart. " >> $logfile
	else
		echo "Status: $exitcode -- $sc $mode $tstart. " >> $logfile
	fi
	
	# Inform the user of error
	if [ ${exitcode} -gt 100 ]; then
		echo "Error: $exitcode -- $sc $mode $tstart."
	fi
done < $list

# Indicate script has finished
echo 'Finished mms_fsm_l3_process_list.sh.'
