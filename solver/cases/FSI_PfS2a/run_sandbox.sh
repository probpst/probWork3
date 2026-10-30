#!/bin/bash
# drive.sh <run dir> [np]: keep an FSI-PfS-2a run going across sandbox restarts.
#  - running already: print status and exit;  HOLD file: do nothing
#  - broke down (stop criterion in the last log): print FAILED and exit (never restarts it)
#  - otherwise resume from the latest restart file (Mode 1) or start (Mode 0); phase 1 = 0.5 ms
#    to step 300 (t = 0.15 s, the release), phase 2 = 1 ms to step 1000 (t = 0.85 s; End_T = steps x Del_T,
#    the solver counts the steps from End_T / Del_T)
d=$1; np=${2:-2}; cd $d || exit 1
[ -f HOLD ] && { echo "HOLD"; exit 0; }
if [ -f .mpirun.pid ] && [ "$(readlink /proc/$(cat .mpirun.pid)/cwd 2>/dev/null)" = "$(pwd -P)" ]; then echo "RUNNING $(tail -1 FSI_PfS2a_monitor.dat 2>/dev/null | awk '{print $1}')"; exit 0; fi
for p in $(pgrep -x probCFD); do
  if [ "$(readlink /proc/$p/cwd)" = "$(pwd -P)" ]; then echo "RUNNING $(tail -1 FSI_PfS2a_monitor.dat | awk '{print $1}')"; exit 0; fi
done
if [ -f log.txt ]; then
  if grep -q 'STOP CRITERION' log.txt; then echo "FAILED $(grep 'CALCULATION STOPPED' log.txt | head -1)"; exit 0; fi
  n=$(ls log_seg*.txt 2>/dev/null | wc -l); mv log.txt log_seg$(printf %02d $((n+1))).txt
fi
last=$(grep -E '^ *[0-9]+ *$' FSI_PfS2a_restart.log 2>/dev/null | tail -1 | tr -d ' ')
endt2=$(cat END 2>/dev/null || echo 1.00)
sl=$(tail -1 FSI_PfS2a_monitor.dat 2>/dev/null | awk '{print $1}')
if [ -n "$sl" ] && awk -v a="$sl" -v b="$endt2" 'BEGIN{exit !(a >= b*1000 - 0.5)}'; then echo "DONE"; exit 0; fi
if [ -z "$last" ]; then mode=0; dt=5.0e-04; endt=0.15
elif [ "$last" -lt 300 ]; then mode=1; dt=5.0e-04; endt=0.15
else mode=1; dt=1.0e-03; endt=$endt2; fi
sed -i "s/^  Mode  .*/  Mode                 = $mode/; s/^ Del_T  .*/ Del_T                 = $dt/; s/^ End_T  .*/ End_T                 = $endt/" probCFD.nml
nohup mpirun --allow-run-as-root -np $np ./probCFD > log.txt 2>&1 &
echo $! > .mpirun.pid
echo "STARTED mode $mode dt $dt End_T $endt from restart step ${last:-none}"
