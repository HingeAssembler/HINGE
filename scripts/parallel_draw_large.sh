#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for i in $(seq 1 1 100)
  do
     echo drawing read $i
     num1=$(ps -ef | grep 'python draw.py' | wc -l)
     num2=$(ps -ef | grep 'LA4Awesome' | wc -l)
     num=$(( $num1 + $num2 ))
     echo $num running
     while [ $num -gt 12 ]
         do 
             sleep 5
             echo waiting
         done
     python draw.py $i &
 done