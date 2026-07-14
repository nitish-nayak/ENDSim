#!/bin/bash

# Make box around text @climagic
function box() { t="$1xxxx";c=${2:-=}; echo ${t//?/$c}; echo "$c $1 $c"; echo ${t//?/$c}; }


MAC="./muon_1gev.mac"

dists=("-4500" "-6500")
angles=("0 -1 0" "-0.5 -0.866 0")

for d in "${dists[@]}"; do
    j=0
    for a in "${angles[@]}"; do
        j=$((j+1))
        NEW_MAC=`echo $MAC | sed -e 's/\.mac//g'`"_d"$d"m_angle"$j".mac"
        cp "$MAC" "$NEW_MAC"
        sed -i "s/pos\/set [0-9-]*/pos\/set $d/g" "$NEW_MAC"
        sed -i "s/vtx\/set mu- [0-9-]* [0-9-]* [0-9-]*/vtx\/set mu- $a/g" "$NEW_MAC"
        sed -i "s/muon_1gev/muon_1gev_d"$d"m_angle"$j"/g" "$NEW_MAC"
        box "running $NEW_MAC"
        ./install/bin/end "$NEW_MAC" > /dev/null 2>&1 &
    done
done
