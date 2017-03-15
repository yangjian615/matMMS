#!/bin/bash
#mms_fg_l2pre.sh
# Usage:
#  For daily srvy output (combines available slow, fast and f128)
# mms_fg_l2pre.sh obs instr srvy YYYYMMDD
#
#  For burst mode processing, there is a one to one correspondence between the input
#  L1A file and the output L2pre file
# mms_fg_l2pre.sh obs instr brst YYYYMMDDhhmmss
#
#  The f128 mode interface is similar to burst if you want 128 S/s output:
# mms_fg_l2pre.sh obs instr f128 YYYYMMDD
#
#  in the above,
#  obs = mms1, mms2, mms3, or mms4
#  instr = afg or dfg
#
# Return status
#  0       Exit with no problems
#  1-100   non-fatal error
#     99       no housekeeping 0x10e for temperatures provided, or data is
#              incomplete
#    100       in nominal operations for srvy, both fast and slow mode L1A inputs
#              are expected, but only one mode found.
#  101-200 Fatal error (no output file produced)
#    126       IDL compilation error
#    127       idl command not found (log created)
#    128       (no log file! msg to stderr) error with mkdir for log file.
#    129       (no log file! msg to stderr) command line usage error
#    130       mms_fg_l2pre IDL function usage error:  improper arguments.
#    139       IDL Segmentation fault!
#    160       reserved for: error reading cal files
#    197       error reading / finding any slow, fast or f128 mode L1A
#              file (for srvy) or
#              error reading /finding specified f128/burst L1A file
#    200       error - Other error (bad cal file, unexpected error, etc.)
#
# see mms_fg_l2pre.pro for list of codes that are explicitly set by that routine.
#
# Output to sdtout
#  Log file name
#  Output file name (or "error" if none produced -- return status 101-200)
#
# Output files created:
#  Log file (unless it is not possible to create a log file!)
#  L2pre file

if [ $# -ne 4 ] && [ $# -ne 5 ];then
  echo 'mms_fg_l2pre: incorrect # of arguments' > /dev/stderr
  echo 'error'
  exit 129
fi

obs=${1}
instr=${2}
mode=${3}
time=${4}
year=${time:0:4}
month=${time:4:2}

date=`date -u +%Y%m%d_%H%M%S`

if [ $# -eq 5 ]; then
	logfile=${5}
else
	logfile="${obs}_${instr}_${mode}_l2pre_${time}_v0.0.0_${date}.log"
fi

logpath=${LOG_PATH_ROOT}/${obs}/${instr}/${mode}/l2pre/${year}/${month}

mkdir -p $logpath
code=$?
if [ $code -ne 0 ] ; then 
  echo 'error'
  exit 128
fi
logpath=${logpath}/${logfile}
echo $logpath

cat > $logpath <<EOF
=========================
Log file for mms_fg_l2pre
Log Filename:     $logpath
Runtime:          $date
Pwd:              `pwd`
Host Information: `uname -a`
LOG_PATH_ROOT:    $LOG_PATH_ROOT
Command Line Arguments:
$0 $@
=========================
IDL_PATH:         $IDL_PATH
IDL_DIR:          $IDL_DIR
calling IDL:
$IDL_DIR/bin/idl mms_fg_l2pre_batch -args $@
=========================
EOF

$IDL_DIR/bin/idl mms_fsm_fgm_l2 -args $@ 2>> $logpath <<EOF
  exit, status=127 ; prevents IDL hanging if the batch file can't be found
EOF

#l2_file="$($IDL_DIR/bin/idl mms_fgm_l2_omb_script -args $@ 2>> $logpath <<EOF
#  exit, status=127 ; prevents IDL hanging if the batch file can't be found
#EOF
#)"

code=$?

# handle case where IDL command is not found:
if [ $code -eq 127 ] ; then echo error; exit 127; fi
exit $code
