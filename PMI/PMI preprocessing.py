# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 10:19:26 2019

@author: cc gong
"""

from collections import defaultdict
import os
import pickle
import tqdm
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
import re
import sys
os.chdir(r'C:\Users\cc gong\Documents\interaction\scripts')


def create_atoms(mol):
    """Create a list of atom (e.g., hydrogen and oxygen) IDs
    considering the aromaticity."""
    atoms = [str(re.findall(pattern='SP\\d{0,1}', string=str(a.GetHybridization()))) for a in mol.GetAtoms()]
    '''
    for a in mol.GetAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], str(a.GetSymbol()))

    for a in mol.GetAromaticAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    '''
    for a in mol.GetAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], a.GetMass())
    atoms = [atom_dict[a] for a in atoms]
    return np.array(atoms)


def create_ijbonddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the tuples of its neighboring node
    and bond (e.g., single and double) IDs."""
    i_jbond_dict = defaultdict(lambda: [])
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond = bond_dict[str(b.GetBondType())]
        i_jbond_dict[i].append((j, bond))
        i_jbond_dict[j].append((i, bond))
    return i_jbond_dict


def creat_hybriddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the atom hybridization ID."""
    hybrid_dict = defaultdict(lambda: len(hybrid_dict))
    hybrids = defaultdict(lambda: [])
    As = [a for a in mol.GetAtoms()]
    for a in As:
        h = hybrid_dict[str(a.GetHybridization())]
        hybrids[a.GetIdx()].append(h)


def extract_fingerprints(atoms, i_jbond_dict, radius):
    """Extract the r-radius subgraphs (i.e., fingerprints)
    from a molecular graph using Weisfeiler-Lehman algorithm."""
    if (len(atoms) == 1) or (radius == 0):
        fingerprints = [fingerprint_dict[a] for a in atoms]
    else:
        nodes = atoms
        i_jedge_dict = i_jbond_dict

        for _ in range(radius):

            """Update each node ID considering its neighboring nodes and edges
            (i.e., r-radius subgraphs or fingerprints)."""
            fingerprints = []
            for i, j_edge in i_jedge_dict.items():
                neighbors = [(nodes[j], edge) for j, edge in j_edge]
                fingerprint = (nodes[i], tuple(sorted(neighbors)))
                fingerprints.append(fingerprint_dict[fingerprint])
            nodes = fingerprints

            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    edge = edge_dict[(both_side, edge)]
                    _i_jedge_dict[i].append((j, edge))
            i_jedge_dict = _i_jedge_dict

    return np.array(fingerprints)


def create_adjacency(mol):
    adjacency = Chem.GetAdjacencyMatrix(mol)
    return np.array(adjacency)


def split_sequence(sequence, ngram):
    for i in range((len(sequence) - ngram + 1)):
        #sequence = '-' + sequence + '='
        words = [word_dict[sequence[i:i + ngram]]]
    return np.array(words)


def dump_dictionary(dictionary, filename):
    with open(filename, 'wb') as f:
        pickle.dump(dict(dictionary), f)


if __name__ == "__main__":
    DATASET, radius, ngram = 'human', 2, 3
    radius, ngram = map(int, [radius, ngram])
    xo = Chem.MolFromPDBFile('./importDat/XO/XO_1fiq.pdb')
    # other=sys.argv[1]
    # other = Chem.MolFromPDBFile('./importDat/XO/'+other+'.pdb')#xo:XO_1fiq;SOD:SOD_1Q0E;CAT:CAT_1DGB;ATP_synthase_5ari;GS_3ng0;GPX:GPX_2ltf;GSHS_2hgs;asp_3wzf;AKT1_3o96;PIK3_4ykn;EGFR_1m17;MAPK_2erk;Uba1_6dc6;UbcH5B_2clw;CDC34_2ob4;HECT_3jvz;UbcH5_CHIP_U-box_complex_2oxq;Structure_of_RING_finger_protein_165_5d0i;PHD_finger_of_human_UHRF1_3zvz
    xoSeq = Chem.MolToSequence(xo)
    #otherSeq = Chem.MolToSequence(other)
    '''
    sys.argv[1]: peptide sequence
    sys.argv[2]: peptide Name
    sys.argv[3]: xo?
    '''
    otherSeq = sys.argv[1]
    mol = {}
    for root, dirs, files in os.walk(r'./importDat/XO'):
        for i in files:
            if i.split(sep='.')[-1] == 'mol' and re.match('\\d+b', i.split(sep='.')[0]) == None:
                path = os.path.join(root, i)
                mol[i.split(sep='.')[0]] = Chem.MolFromMolFile(path)
            else:
                pass
    mol['19a'] = Chem.MolFromMolFile('./importDat/XO/19a.mol')
    #mol['19b'] = Chem.MolFromMolFile('./importDat/XO/19b.mol')
    data_list_xo = {}
    for k, v in mol.items():
        if re.match('\\d+b', k) is not None:
            i = 0
        elif k in [str(n) + 'a' for n in range(1, 17)]:
            i = 1
        elif k in [str(n) + 'a' for n in range(17, 21)]:
            i = 0
        else:
            i = 1
        data_list_xo[k] = Chem.MolToSmiles(v) + ' ' + xoSeq + ' ' + str(i)
    if sys.argv[3] == 'xo':
        data_list = data_list_xo
    else:
        data_list_other = {}
        for k, v in mol.items():
            i = 0
            data_list_other[k] = Chem.MolToSmiles(v) + ' ' + otherSeq + ' ' + str(i)
        data_list = data_list_other

    N = len(data_list)
    atom_dict = defaultdict(lambda: len(atom_dict))
    bond_dict = defaultdict(lambda: len(bond_dict))
    fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
    edge_dict = defaultdict(lambda: len(edge_dict))
    word_dict = defaultdict(lambda: len(word_dict))
    Smiles, compounds, adjacencies, proteins, interactions, ID = '', [], [], [], [], []
    for index, data in tqdm.tqdm(data_list.items()):
        smiles, sequence, interaction = data.strip().split()
        Smiles += smiles + '\n'
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))  # Consider hydrogens.
        atoms = create_atoms(mol)
        i_jbond_dict = create_ijbonddict(mol)
        fingerprints = extract_fingerprints(atoms, i_jbond_dict, radius)
        ID.append(index)
        compounds.append(fingerprints)
        adjacency = create_adjacency(mol)
        adjacencies.append(adjacency)
        words = split_sequence(sequence, ngram)
        proteins.append(words)
        interactions.append(np.array([float(interaction)]))
    dir_input = (r'./validation/inputDat/' + 'radius_' + str(radius) + '_ngram_' + str(ngram) + sys.argv[2] + '/')
    os.makedirs(dir_input, exist_ok=True)
    with open(dir_input + 'Smiles.txt', 'w') as f:
        f.write(Smiles)
    np.save(dir_input + 'index', ID)
    np.save(dir_input + 'compounds', compounds)
    np.save(dir_input + 'adjacencies', adjacencies)
    np.save(dir_input + 'proteins', proteins)
    np.save(dir_input + 'interactions', interactions)
    dump_dictionary(fingerprint_dict, dir_input + 'fingerprint_dict.pickle')
    dump_dictionary(word_dict, dir_input + 'word_dict.pickle')
    print('\nThe preprocess of ' + DATASET + ' dataset has finished!')
