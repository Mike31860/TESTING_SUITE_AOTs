#!/bin/sh
GIT_DIR=$PWD
F_LOG=$GIT_DIR/3mm.txt
printf "try" > $F_LOG
./serial_3mm | grep 'Time in seconds'>> $F_LOG
./parallel_3mm | grep 'Time in seconds'>> $F_LOG
echo "<----------------------------->"
