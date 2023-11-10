#!/bin/bash 
rm=5.
mg=(1 10 30 50 100 )
for ((i = 0; i < ${#mg[@]}; i++))
do
    #a="scale=2; 5.0 + 0.1 * $i"|bc
    #z="scale=2; $rm + ${mg[$i]}"|bc
    g++ -DTAU=${mg[$i]} -O3 rststwkannp_mr.cpp -o ana$i.exe
    ./ana$i.exe
    echo "${mg[$i]} $i"
done