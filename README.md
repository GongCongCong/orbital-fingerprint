# orbital-fingerprint

 We used molecular characterization (hybrid orbitals and electron distribution) to construct a mathematical model to evaluate the ability to form hydrogen bonds and named it with lone-pair electron (LPE).To comprehensively explore molecular interaction, using molecular properties including the hybrid orbital forms of each atom, atomic bonding types, and the adjacency list of atoms, we constructed a molecular fingerprint and named it orbital fingerprint. Then a graph neural network (GNN) for compounds was constructed using orbital fingerprint. A prediction for molecular interaction (PMI) was performed by combining this GNN and a convolutional neural network (CNN) for proteins.

## LPE

![LPE formula](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;\huge&space;LPE=&space;\frac{\sum^{n}_{1}N_i^e}{N_{sp^2}&space;&plus;&space;N_{{sp}^3}})

<div align=center><img src="https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;\huge&space;LPE=&space;\frac{\sum^{n}_{1}N_i^e}{N_{sp^2}&space;&plus;&space;N_{{sp}^3}}" width="150" height="150" alt="LPE formula" /></div>

where 
![N_i^e](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;N_i^e) is the number of lone paired electron of the 
![i](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;i)-th atom; 
![N_{sp^2}](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;N_{sp^2}) and 
![N_{sp^3}](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;N_{sp^3}) are the numbers of 
![sp^2](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;sp^2)- and 
![sp^3](https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;sp^3)-hybridized atoms, respectively. 

## PMI

Here, we reference the CPI approach (https://github.com/masashitsubaki/CPI_prediction) to predict compound-protein interaction. But the input information of PMI for the compound is the orbital fingerprint.
