#!/bin/sh
GIT_DIR=$PWD
F_LOG=$GIT_DIR/atax.txt
printf "try" > $F_LOG
./serial_atax | grep 'Time in seconds'>> $F_LOG
./parallel_atax | grep 'Time in seconds'>> $F_LOG
echo "<----------------------------->"
