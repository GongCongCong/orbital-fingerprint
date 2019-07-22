# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 22:05:08 2019

@author: cc gong
"""
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import os
import pickle
import numpy as np
import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
os.chdir(r'C:\Users\cc gong\Documents\interaction\scripts')
setting = r'human--radius2--ngram3--dim10--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration100'
modelPath = r'./exportDat/outputDat/model/' + setting + '.model'


DATASET = 'human'
dim = 10
radius = 2
ngram = 3
layer_gnn = 3
layer_cnn = 3
window = 5
layer_output = 3
lr = 1e-4
lr_decay = 0.5
weight_decay = -6.
decay_interval = 10
iteration = 100

if torch.cuda.is_available():
    device = torch.device('cuda')
    print('The code uses GPU...')
else:
    device = torch.device('cpu')
    print('The code uses CPU!!!')


def load_tensor(file_name, dtype):
    return [dtype(d).to(device) for d in np.load(file_name + '.npy')]


def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)


dir_input = (r'./validation/inputDat/' + 'radius_' + str(radius) + '_ngram_' + str(ngram) + sys.argv[1] + '/')
compounds = load_tensor(dir_input + 'compounds', torch.LongTensor)
adjacencies = load_tensor(dir_input + 'adjacencies', torch.FloatTensor)
proteins = load_tensor(dir_input + 'proteins', torch.LongTensor)
interactions = load_tensor(dir_input + 'interactions', torch.LongTensor)
fingerprint_dict = load_pickle(dir_input + 'fingerprint_dict.pickle')
word_dict = load_pickle(dir_input + 'word_dict.pickle')
ID = np.load(dir_input + 'index.npy')
ID = ID.tolist()
n_fingerprint = len(fingerprint_dict)
n_word = len(word_dict)


class CPIPre(nn.Module):
    def __init__(self):
        super(CPIPre, self).__init__()
        self.embed_fingerprint = nn.Embedding(n_fingerprint, dim)
        self.embed_word = nn.Embedding(n_word, dim)
        self.W_gnn = nn.ModuleList([nn.Linear(dim, dim) for _ in range(layer_gnn)])
        self.W_cnn = nn.ModuleList([nn.Conv2d(
            in_channels=1, out_channels=1, kernel_size=2 * window + 1,
            stride=1, padding=window) for _ in range(layer_cnn)])
        self.W_attention = nn.Linear(dim, dim)
        self.W_out = nn.ModuleList([nn.Linear(2 * dim, 2 * dim)
                                    for _ in range(layer_output)])
        self.W_interaction = nn.Linear(2 * dim, 2)

    def gnn(self, xs, A, layer):
        for i in range(layer):
            hs = torch.relu(self.W_gnn[i](xs))
            xs = xs + torch.matmul(A, hs)
        # return torch.unsqueeze(torch.sum(xs, 0), 0)
        return torch.unsqueeze(torch.mean(xs, 0), 0)

    def attention_cnn(self, x, xs, layer):
        """The attention mechanism is applied to the last layer of CNN."""

        xs = torch.unsqueeze(torch.unsqueeze(xs, 0), 0)
        for i in range(layer):
            xs = torch.relu(self.W_cnn[i](xs))
        xs = torch.squeeze(torch.squeeze(xs, 0), 0)

        h = torch.relu(self.W_attention(x))
        hs = torch.relu(self.W_attention(xs))
        weights = torch.tanh(F.linear(h, hs))
        ys = torch.t(weights) * hs

        # return torch.unsqueeze(torch.sum(ys, 0), 0)
        return torch.unsqueeze(torch.mean(ys, 0), 0)

    def forward(self, inputs):
        fingerprints, adjacency, words = inputs

        """Compound vector with GNN."""
        fingerprint_vectors = self.embed_fingerprint(fingerprints)
        compound_vector = self.gnn(fingerprint_vectors, adjacency, layer_gnn)

        """Protein vector with attention-CNN."""
        word_vectors = self.embed_word(words)
        protein_vector = self.attention_cnn(compound_vector,
                                            word_vectors, layer_cnn)

        """Concatenate the above two vectors and output the interaction."""
        cat_vector = torch.cat((compound_vector, protein_vector), 1)
        for j in range(layer_output):
            cat_vector = torch.relu(self.W_out[j](cat_vector))
        interaction = self.W_interaction(cat_vector)

        return interaction

    def __call__(self, data, train=True):

        inputs, correct_interaction = data[:3], data[-1]
        predicted_interaction = self.forward(inputs)

        if train:
            loss = F.cross_entropy(predicted_interaction, correct_interaction)
            return loss
        else:
            correct_labels = correct_interaction.to('cpu').data.numpy()
            ys = F.softmax(predicted_interaction, 1).to('cpu').data.numpy()
            predicted_labels = list(map(lambda x: np.argmax(x), ys))
            predicted_scores = list(map(lambda x: x[1], ys))
            return correct_labels, predicted_labels, predicted_scores


def shuffle_dataset(dataset, seed):
    np.random.seed(seed)
    np.random.shuffle(dataset)
    return dataset


model = torch.load(modelPath)
dataset = list(zip(compounds, adjacencies, proteins, ID, interactions))
dataset = shuffle_dataset(dataset, 1234)

res = {}

for data in dataset:
    res[data[3]] = [model(data, train=False)]

file_pre = r'./exportDat/outputDat/result/prediction-' + sys.argv[1] + '-' + setting + '.txt'
with open(file_pre, 'w') as f:
    f.write('Name'.center(10, ' ') + '\tCorrect\tPredict\tTorF' + '\n')

print('Name'.center(10, ' '), '\t' + 'Correct\t' + 'Predict\t')
true = 0
false = []
for k, v in res.items():
    n = k.center(10, ' ')
    c = v[0][0].tolist()
    p = v[0][1]
    t = ['T' if c == p else 'F']
    if t[0] == 'T':
        true += 1
    if t[0] == 'F':
        false.append(k)
    with open(file_pre, 'a') as f:
        f.write(str(n) + '\t' + str(c[0]) + '\t' + str(p[0]) + '\t' + str(t[0]) + '\n')
    print('{}\t{}\t{}\t{}'.format(n, c, p, t))
true / len(res)
