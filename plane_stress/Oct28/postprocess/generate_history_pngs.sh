#!/bin/bash

SNOPTOUTS=*SNOPT*.out
IPOPTOUTS=*IPOPT*.out
PAROPTOUTS=*ParOpt*.tr

for f in $SNOPTOUTS
  do 
    echo "Generating history plot for $f ..."
    snopt_plot $f --save
done

for f in $IPOPTOUTS
do
  echo "Generating history plot for $f ..."
  ipopt_plot $f --save
done

for f in $PAROPTOUTS
do
  echo "Generating history plot for $f ..."
  paropt_plot $f --save
done

echo "Done generating all history plots!"
