#!/bin/bash 
mg=(5.5 5.6 5.7 5.8 5.9 5.5 )
for ((i = 0; i < ${#mg[@]}; i++))
do
    icc -std=c++14 -DR=${mg[$i]} -O3 rststwkannp_mr.cpp -o ana$i.out
    qsub twowstw$i.sh
    echo "icc -std=c++14 -DRbit=${mg[$i]} -O3 rststwkan.cpp -o ana$i.out"
done