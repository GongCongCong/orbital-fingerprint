# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 10:52:14 2019

@author: cc gong
"""

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd
import os
from rdkit.Chem import Descriptors
import re
# In[class structrue]:
class Structrue():
    def __init__(self, name = None, smi = None, mol = None):
        if not mol == None:
            mol = mol
        elif not smi == None:
            mol = Chem.MolFromSmiles(smi)
        else:
            raise NameError('mol or smi missed'.center(40,'='))
        self.mol = mol
    def cal_info(self):
        symbol = []
        Implicit = []
        Explicit = []
        total = []
        AtomicNum = []
        hybrid = []
        e = []
        for atom in self.mol.GetAtoms():
            symbol.append(atom.GetSymbol())
            Implicit.append(atom.GetImplicitValence())
            Explicit.append(atom.GetExplicitValence())
            total.append(atom.GetTotalValence())
            AtomicNum.append(atom.GetAtomicNum())
            hybrid.append(str(atom.GetHybridization()))
            if atom.GetSymbol() in ['O','S']:
                e.append(2)
            elif atom.GetSymbol() in ['N','P']:
                e.append(3)
            else:
                e.append(0)
        dict={'symbol':symbol,'Implicit':Implicit,'Explicit':Explicit,'total':total,'AtomicNum':AtomicNum,'hybridization':hybrid,'e':e}
        print(r'completed!'.center(40,'+'))
        return pd.DataFrame(dict)
    def LEA(self):
        """
        Nsp2: numbers of sp2-hybridized atoms
        Nsp3: numbers of sp3-hybridized atoms
        Ne: Lone pair electron number
        Na: Lone pair electron conjugate number
        Nb: Chemical bond number
        Nc: Number of conjugated bonds
        N: atom number
        """
        Nsp2 = 0
        Nsp3 = 0
        Ne = 0
        Na = 0
        Nb = 0
        conjugated = []
        N = self.mol.GetNumAtoms()
        Nb = self.mol.GetNumBonds()
        for atom in self.mol.GetAtoms():
            bond = atom.GetBonds()
            for i in bond:
                conjugated.append(i.GetIsConjugated())
            if str(atom.GetHybridization()) == 'SP2':
                Nsp2 += 1
            elif str(atom.GetHybridization()) == 'SP3':
                Nsp3 += 1
            else:
                continue
            if atom.GetSymbol() in ['O','S']:
                Ne += 2
                for i in atom.GetNeighbors():
                    temp = []
                    for j in i.GetBonds():
                        temp.append(str(j.GetBondType()))
                    if temp ==['AROMATIC', 'AROMATIC', 'SINGLE']:
                        Na += 1
            elif atom.GetSymbol() in ['N','P']:
                for i in atom.GetNeighbors():
                    temp = []
                    for j in i.GetBonds():
                        temp.append(str(j.GetBondType()))
                    if temp ==['AROMATIC', 'AROMATIC', 'SINGLE']:
                        Na += 1
                Ne += 3
            else:
                Ne = Ne
        Nc = sum(conjugated)
        LEA = Ne/(Nsp3+Nsp2)
        return LEA
def LEAFromSeq(seq):
    seq = re.sub('[^A-Z]', '', seq)
    struc = Structrue(mol = Chem.MolFromSequence(seq,flavor=0))
    lea = struc.LEA()
    return lea
# In[读取20个氨基酸的结构]
aa = pd.read_table(r'C:\Users\cc gong\Documents\interaction\scripts\importDat\buildinDat\amino acid.txt',sep=',').loc[:,'ID_mono':'亲疏水']
path = r'C:\Users\cc gong\Documents\interaction\scripts\buildinDat'
aaMolPeps = {}
for root, dirs, files in os.walk(path, topdown=False):
    for name in files:
        if str(name).rsplit(sep='.')[1] == 'sdf':
            mol = Chem.MolFromMolFile(os.path.join(root, name))
            aaMolPeps[name.split(sep='.')[0]] = mol
# In[计算20氨基酸]
aaWLEA = []
for k in aa.ID_mono:
    aaWLEA.append(Structrue(mol=aaMolPeps[k]).LEA())
pd.DataFrame({'name':aa.name,'w':aa.loc[:,'亲疏水'],'LEA':aaWLEA})
# In[plot a,f,i,l,w,y]
for a in ['A','F','I','L','W','Y']:
    print(a)
    display(aaMolPeps[a])
# In[读取peptide结构]
path = r'C:\Users\cc gong\Documents\interaction\scripts\peptide lps'
molPeps = {}
for root, dirs, files in os.walk(path, topdown=False):
    for name in files:
        if str(name).rsplit(sep='.')[1] == 'sdf':
            mol = Chem.MolFromMolFile(os.path.join(root, name))
            molPeps[name] = mol
# In[计算]
pepStructure = {}
pepLEA = {}
for k,v in molPeps.items():
    pepStructure[k] = Structrue(name = k, mol=v)
for k,v in pepStructure.items():
    pepLEA[k] = v.LEA()
pepLEA
# In[BIOGRID]
path = r'C:\Users\cc gong\Documents\interaction\scripts\peptide lps\BIOGRID-ALL-3.5.171-20190411'
bioGrid = pd.read_table(os.path.join(path,'BIOGRID-ALL-3.5.171.tab.txt'), sep='\t', comment ='#',nrows=1000)
bioGrid = bioGrid.loc[:,['OFFICIAL_SYMBOL_A','OFFICIAL_SYMBOL_B','EXPERIMENTAL_SYSTEM']]
Sequence = pd.read_table(os.path.join(path,'Sequence.txt'), sep='\t', comment ='#')
Sequence = Sequence.loc[:, ['peptide','hgnc_symbol']]
# In[]
mol_A ={}
mol_B ={}
for i in range(20,40):
    print('开始第{}行'.format(int(i)+1).center(40,'#'))
    A = bioGrid.loc[i,'OFFICIAL_SYMBOL_A']
    B = bioGrid.loc[i,'OFFICIAL_SYMBOL_B']
    print('{}有{}个蛋白序列！\n'.format(A,len(Sequence.loc[Sequence.hgnc_symbol==A,'peptide'])))
    print('{}有{}个蛋白序列！\n'.format(B,len(Sequence.loc[Sequence.hgnc_symbol==B,'peptide'])))
    AIdx = 1
    BIdx = 1
    for As in Sequence.loc[Sequence.hgnc_symbol==A,'peptide']:
        try:
            s = re.sub('[^A-Z]+', '', As)
            print('计算{}中的第{}个蛋白序列：\n{}\n'.format(A, AIdx,s))
            lea = LEAFromSeq(s)
            mol_A[A+'_'+str(AIdx)] = lea
            AIdx +=1
        except NameError:
            print('计算{}中的第{}个蛋白序列出错：\n{}\n'.format(A, AIdx,s))
            pass
        continue
    for Bs in Sequence.loc[Sequence.hgnc_symbol==B,'peptide']:
        try:
            s = re.sub('[^A-Za-z]+', '', Bs)
            print('计算{}中的第{}个蛋白序列：\n{}\n'.format(B, BIdx,s))
            lea = LEAFromSeq(s)
            mol_B[B+'_'+str(BIdx)] = lea
            BIdx +=1
        except NameError:
            print('计算{}中的第{}个蛋白序列出错：\n{}\n'.format(B, BIdx,s))
            pass
        continue

mol_A_1 = pd.DataFrame.from_dict(mol_A,orient='index').rename(columns={0:'lea'})
mol_B_1 = pd.DataFrame.from_dict(mol_B,orient='index').rename(columns={0:'lea'})
mol_A_1.to_excel('C:\\Users\\cc gong\\Documents\\interaction\\scripts\\exportDat\\A_lea_from_bioGrid.xlsx', index_label='symbol', sheet_name='A')
mol_B_1.to_excel('C:\\Users\\cc gong\\Documents\\interaction\\scripts\\exportDat\\B_lea_from_bioGrid.xlsx', index_label='symbol', sheet_name='B')
#bioGrid.to_excel('C:\\Users\\cc gong\\Documents\\interaction\\scripts\\exportDat\\bioGrid.xlsx', index=False)
# In[]
from tqdm import tqdm
app = Sequence.loc[Sequence.hgnc_symbol=='APP', ['peptide']].values
lea = {}
mol = {}
for i in tqdm(range(len(app))):
  try:
    seq = re.sub('[^A-Z]', '', app[i].tolist()[0])
    mol[i] = Chem.MolFromSequence(seq,flavor=0)
    struc = Structrue(mol = mol)
    l = struc.LEA()
    lea[i] = [l,seq]
  except NameError:
    pass
