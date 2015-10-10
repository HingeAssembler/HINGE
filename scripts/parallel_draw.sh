#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for i in $(seq 4000 1 20000)
  do
     echo drawing read $i
     num1=$(ps -ef | grep 'python draw.py' | wc -l)
     num2=$(ps -ef | grep 'LA4Awesome' | wc -l)
     num=$(( $num1 + $num2 ))
     echo $num running
     while [ $num -gt 60 ]
         do 
             sleep 5
             echo waiting, $num running
             num1=$(ps -ef | grep 'python draw.py' | wc -l)
             num2=$(ps -ef | grep 'LA4Awesome' | wc -l)
             num=$(( $num1 + $num2 ))             
         done
     python draw2.py $i &
 done
