#!/bin/sh
GIT_DIR=$PWD
F_LOG=$GIT_DIR/2mm.txt
printf "try" > $F_LOG
./serial_2mm | grep 'Time in seconds'>> $F_LOG
echo "<----------------------------->"
