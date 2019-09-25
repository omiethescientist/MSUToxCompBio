#!/usr/bin/env python
# coding: utf-8

# In[6]:


#Credit to Patrick Ng for his DNA2Vec Package
import twobitreader as tbr
import re 
import numpy as np


# In[7]:


def chr_strings(genome):
    """
    Create a string generator for chromosome 
    """
    for chrom in genome:
        yield genome[chrom][:]

def rev_comp(seq):
    """
    Get the DNA reverse complements
    """
    comps = {"a":"t", "t":"a", "g":"c", "c":"g"}
    seq = seq.lower()
    return "".join([comps[i] for i in seq[::-1]])

def get_pseudo_tads(string, motif):
    """
    split the DNA according to some motif
    The motif is written in the form of a regex 
    """
    x = re.split(motif, string.lower())
    x = filter(bool, x)
    for i in x:
        if "n" not in i:
            if np.random.rand() < 0.05:
                yield rev_comp(i)
            else:
                yield(i)
            
def kmerFragmenter(seq, k):
    return " ".join([seq[i:i+k] for i in range(len(seq)-k)])


# In[8]:


genome = tbr.TwoBitFile("../data/hg38.2bit")
motif = r"agggggc"
k = 8
file = "../data/humanGenomeNLP/"
n = 0
chroms = chr_strings(genome)
for i in chroms:
    for j in get_pseudo_tads(str(i), motif):
        with open(file + "test_" + str(n), "w+") as f:
            n += 1
            f.write(kmerFragmenter(j, k))


# In[ ]:




