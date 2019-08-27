# orbital-fingerprint

 We used molecular characterization (hybrid orbitals and electron distribution) to construct a mathematical model to evaluate the ability to form hydrogen bonds and named it with lone-pair electron (LPE).To comprehensively explore molecular interaction, using molecular properties including the hybrid orbital forms of each atom, atomic bonding types, and the adjacency list of atoms, we constructed a molecular fingerprint and named it orbital fingerprint. Then a graph neural network (GNN) for compounds was constructed using orbital fingerprint. A prediction for molecular interaction (PMI) was performed by combining this GNN and a convolutional neural network (CNN) for proteins.

## LPE

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

## The calculation of the lone-pair electron (LPE) index
The molecule structure was converted to simplified molecular-input line-entry system (SMILES) format using RDKit (Open-source cheminformatics; http://www.rdkit.org). And we used the Chem module of RDKit to calculate the numbers of lone pair electrons and the numbers of the sp2- and sp3-orbitals of each molecules. Then the LPE index was calculated using the formula: 
<div align=center><img src="https://latex.codecogs.com/svg.latex?\inline&space;\bg_white&space;\fn_cs&space;\huge&space;LPE=&space;\frac{\sum^{n}_{1}N_i^e}{N_{sp^2}&space;&plus;&space;N_{{sp}^3}}" width="150" height="150" alt="LPE formula" /></div>

## The construction of the orbital fingerprint
The molecule graph is represented as g=(V,E) to the molecular graph, where v_i∈V is the hybridization type (e.g., sp2 and sp3) of i-th atom and ε_ij∈E is the chemical bond type (e.g., single and double) between the i-th and j-th atoms(22) (Figure 2a & Figure S2a). Since biomolecules are kind of complex, the above graph is far enough to decipher tiny differentiation between molecules owning close structure, the r-radius subgraphs as an effective make-up strategy was applied as well(15). The r value is defined as the number of hops from a core atom. In the present study, r = 2 was used the number of time steps in a graph neural network(23) to predict the molecule interaction because it has already proved that the depth of two steps. To be more precise, we described the hybridization types of neighboring atoms and chemical bonds types within the radius of the two atoms from the i-th vertex as N_v and N_e, respectively. Then we used the different non-negative integers to represent the different types of N_v and N_e, and we assigned an embedding to characterize a specific molecule. For example, [CH2=CH-CH2-CH3] could be described by different r-radius subgraphs (r = 2) as “sp2=sp2-sp3”, “sp2=sp2-sp3-sp3”, “sp2=sp2-sp3-sp3”, and “sp2-sp3-sp3”. Based on this, we dispatched four non-negative integers 0, 1, 1, and 2 to the above four orbital expressions, respectively. Therefore, we could use the vector (0, 1, 1, 2) to represent this molecule (Figure S2a). To be more precisely characterize the molecule, the adjacency list of atoms was also considered in the molecular graph. Therefore, by graphing each atomic hybrid orbital, we define a very different representation of molecule, which we name it as an orbital fingerprint. The detailed construction steps of the orbital fingerprint are shown in https://github.com/GongCongCong/orbital-fingerprint. 
