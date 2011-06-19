#!/bin/bash

OS=`uname -s`

echo "Exponential Decay"
./curvefit -f examples/expdecay_ex.txt -m expdecay

echo "Gaussian"
./curvefit -f examples/gaussian_ex.txt -m gaussian -c 40 -d 5

echo "Hill equation"
./curvefit -f examples/hill_ex.txt -m hill

echo "Dose response (ic50)"
./curvefit -f examples/ic50_ex.txt -m ic50

echo "Michaelis-Mentin"
./curvefit -f examples/mm_ex.txt -m mm

echo "Modified sine wave"
./curvefit -f examples/modsin_ex.txt -m modsin -c 2.0

for file in examples/*_script.plt
do
	gnuplot < $file
done

if [ "$OS" == "Darwin" ]
then
	`open examples/*.png`
elif [ "$OS" == "Linux" ]
then
	`gnome-open examples/*.png`
fi