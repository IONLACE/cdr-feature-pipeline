#!/bin/bash

for i  in 3rjq 4c59 4krl 4wen 5c1m 5f21 5vxl 6fuz 6n48 6obm 6rqm 6x07 6yu8 7a6o 7ar0 7bts 7d4b 7rnn 7th3 7w1s 7wki 8b7w 8gjr 8pij 8qf4 8s61 8uoa 8uw9 8xps 8y9t 8z8v 9fww 9fzd 9goe 9ho4 9j3l/
do

 dockq=$(grep iRMSD $i/idx0Repair.dockq | awk '{print $17"," $2"," $4"," $6"," $8}') 
 dG=$(grep "Total    " $i/0_binding.log | awk '{print $NF}')

echo $dockq,$dG

 dockq=$(grep iRMSD $i/idx1Repair.dockq | awk '{print $17"," $2"," $4"," $6"," $8}') 
 dG=$(grep "Total    " $i/1_binding.log | awk '{print $NF}')

echo $dockq,$dG

 dockq=$(grep iRMSD $i/idx2Repair.dockq | awk '{print $17"," $2"," $4"," $6"," $8}') 
 dG=$(grep "Total    " $i/2_binding.log | awk '{print $NF}')

echo $dockq,$dG

 dockq=$(grep iRMSD $i/idx3Repair.dockq | awk '{print $17"," $2"," $4"," $6"," $8}') 
 dG=$(grep "Total    " $i/3_binding.log | awk '{print $NF}')

echo $dockq,$dG

 dockq=$(grep iRMSD $i/idx4Repair.dockq | awk '{print $17"," $2"," $4"," $6"," $8}') 
 dG=$(grep "Total    " $i/4_binding.log | awk '{print $NF}')

echo $dockq,$dG

# For crystal structure, we set the dockq score to 1 if iRMSD is 0, else 0. We set the other scores to 0 since they are not defined for the crystal structure.
# We also set the binding energy to the value from the binding.log file for the crystal structure.
 dockq_crys=$(grep iRMSD $i/idx0Repair.dockq | awk '{print $21", 1, 0, 0, 0"}')
  dG=$(grep "Total    " $i/binding.log | awk '{print $NF}')

echo  $dockq_crys,$dG

done
