#!/usr/bin/env python
# coding: utf-8

# # Working with DNA2Vec to create Sequence Embeddings
# ## Ng et. al 2017
# We will first test this approach with the DNA2Vec to see if it has any merit.
# In[49]:


### Attic ####

## Utils ##
import random
import string
import resource
import logbook
import arrow
import numpy as np
import os

def split_Xy(df, y_colname='label'):
    X = df.drop([y_colname], axis=1)
    y = df[y_colname]
    return (X, y)

def shuffle_tuple(tup, rng):
    lst = list(tup)
    rng.shuffle(lst)
    return tuple(lst)

def shuffle_dataframe(df, rng):
    """
    this does NOT do in-place shuffling
    """
    return df.reindex(rng.permutation(df.index))

def random_str(N):
    return ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for _ in range(N))

def memory_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6

def estimate_bytes(filenames):
    return sum([os.stat(f).st_size for f in filenames])

def get_output_fileroot(dirpath, name, postfix):
    return '{}/{}-{}-{}-{}'.format(
        dirpath,
        name,
        arrow.utcnow().format('YYYYMMDD-HHmm'),
        postfix,
random_str(3))

## Time Benchmark ##
import time
import logbook

class Benchmark():
    def __init__(self):
        self.time_wall_start = time.time()
        self.time_cpu_start = time.process_time()
        self.logger = logbook.Logger(self.__class__.__name__)

    def diff_time_wall_secs(self):
        return (time.time() - self.time_wall_start)

    def print_time(self, label=''):
        self.logger.info("%s wall=%.3fm cpu=%.3fm" % (
            label,
            self.diff_time_wall_secs() / 60.0,
            (time.clock() - self.time_cpu_start) / 60.0,
        ))

## Tee ##
import sys

class Tee(object):
    def __init__(self, fptr):
        self.file = fptr

    def __enter__(self):
        self.stdout = sys.stdout
        sys.stdout = self

    def __exit__(self, exception_type, exception_value, traceback):
        sys.stdout = self.stdout

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

## Bio Util ##
import numpy as np
import os
from Bio.Seq import Seq

NUCLEOTIDES = 'ACGT'

class Tuple4:
    def __init__(self, pos1, pos2, neg1, neg2):
        self.pos1 = pos1
        self.pos2 = pos2
        self.neg1 = neg1
        self.neg2 = neg2

def determine_out_filename(output_dir, fileroot, mode, extension='txt'):
    return os.path.join(output_dir, '{}.{}.{}'.format(fileroot, mode, extension))

def create_tuple4(kmer1, kmer2, kmer1_neg, kmer2_neg):
    """
    all inputs are list of single nucleotides, e.g. ['A', 'A', 'C']
    """
    return Tuple4(
        ''.join(kmer1),
        ''.join(kmer2),
        ''.join(kmer1_neg),
        ''.join(kmer2_neg))

def insert_snippet(seq, snippet, idx):
    """
    idx: 0 <= idx <= len(seq)]
    """
    split1 = seq[:idx]
    split2 = seq[idx:]
    return split1 + snippet + split2

def pairwise_key(v1, v2):
    return '{}:{}'.format(v1, v2)

def rand_kmer(rng, k_low, k_high=None):
    """
    k_low and k_high are inclusive
    """
    if k_high is None:
        k_high = k_low
    k_len = rng.randint(k_low, k_high + 1)
    return ''.join([NUCLEOTIDES[x] for x in rng.randint(4, size=k_len)])

def rand_nt(rng):
    return NUCLEOTIDES[rng.randint(4)]

def generate_revcompl_pair(k_low, k_high=None, rng=None):
    # TODO make params k_high and rng be required
    if k_high is None:
        k_high = k_low
    if rng is None:
        rng = np.random
    kmer = rand_kmer(rng, k_low, k_high)
    return (kmer, revcompl(kmer))

def revcompl(kmer):
    return str(Seq(kmer).reverse_complement())

def generate_1nt_mutation_4tuple(rng, k_len):
    kmer1 = list(rand_kmer(rng, k_len, k_len))
    kmer2 = list(rand_kmer(rng, k_len, k_len))

    idx = rng.randint(len(kmer1))
    original_nt = kmer1[idx]
    mutate_nt = rand_nt(rng)

    kmer1_neg = list(kmer1)
    kmer1_neg[idx] = mutate_nt

    kmer2[idx] = mutate_nt
    kmer2_neg = list(kmer2)
    kmer2_neg[idx] = original_nt

    return create_tuple4(kmer1, kmer2, kmer1_neg, kmer2_neg)


# In[50]:


### Generators ####
import logbook
import re
from Bio import SeqIO
from itertools import islice
import numpy as np

def remove_empty(str_list):
    return filter(bool, str_list)  # fastest way to remove empty string

class SeqFragmenter:
    """
    Split a sequence into small sequences based on some criteria, e.g. 'N' characters
    """
    def __init__(self):
        pass

    def get_acgt_seqs(self, seq):
        return remove_empty(re.split(r'[^ACGTacgt]+', str(seq)))

class SlidingKmerFragmenter:
    """
    Slide only a single nucleotide
    """
    def __init__(self, k_low, k_high):
        self.k_low = k_low
        self.k_high = k_high

    def apply(self, rng, seq):
        return [seq[i: i + rng.randint(self.k_low, self.k_high + 1)] for i in range(len(seq) - self.k_high + 1)]

class DisjointKmerFragmenter:
    """
    Split a sequence into kmers
    """
    def __init__(self, k_low, k_high):
        self.k_low = k_low
        self.k_high = k_high

    @staticmethod
    def random_chunks(rng, li, min_chunk, max_chunk):
        """
        Both min_chunk and max_chunk are inclusive
        """
        it = iter(li)
        while True:
            head_it = islice(it, rng.randint(min_chunk, max_chunk + 1))
            nxt = '' . join(head_it)

            # throw out chunks that are not within the kmer range
            if len(nxt) >= min_chunk:
                yield nxt
            else:
                break

    def apply(self, rng, seq):
        seq = seq[rng.randint(self.k_low):]  # randomly offset the beginning to create more variations
        return list(DisjointKmerFragmenter.random_chunks(rng, seq, self.k_low, self.k_high))

class SeqMapper:
    def __init__(self, use_revcomp=True):
        self.use_revcomp = use_revcomp

    def apply(self, rng, seq):
        seq = seq.upper()
        if self.use_revcomp and rng.rand() < 0.5:
            return seq.reverse_complement()
        else:
            return seq

class SeqGenerator:
    def __init__(self, filenames, nb_epochs, seqlen_ulim=5000):
        self.filenames = filenames
        self.nb_epochs = nb_epochs
        self.seqlen_ulim = seqlen_ulim
        self.logger = logbook.Logger(self.__class__.__name__)
        self.logger.info('Number of epochs: {}'.format(nb_epochs))

    def filehandle_generator(self):
        for curr_epoch in range(self.nb_epochs):
            for filename in self.filenames:
                with open(filename) as file:
                    self.logger.info('Opened file: {}'.format(filename))
                    self.logger.info('Memory usage: {} MB'.format(memory_usage()))
                    self.logger.info('Current epoch: {} / {}'.format(curr_epoch + 1, self.nb_epochs))
                    yield file

    def generator(self, rng):
        for fh in self.filehandle_generator():
            # SeqIO takes twice as much memory than even simple fh.readlines()
            for seq_record in SeqIO.parse(fh, "fasta"):
                whole_seq = seq_record.seq
                self.logger.info('Whole fasta seqlen: {}'.format(len(whole_seq)))
                curr_left = 0
                while curr_left < len(whole_seq):
                    seqlen = rng.randint(self.seqlen_ulim // 2, self.seqlen_ulim)
                    segment = seq_record.seq[curr_left: seqlen + curr_left]
                    curr_left += seqlen
                    self.logger.debug('input seq len: {}'.format(len(segment)))
                    yield segment

class KmerSeqIterable:
    def __init__(self, rand_seed, seq_generator, mapper, seq_fragmenter, kmer_fragmenter, histogram):
        self.logger = logbook.Logger(self.__class__.__name__)
        self.seq_generator = seq_generator
        self.mapper = mapper
        self.kmer_fragmenter = kmer_fragmenter
        self.seq_fragmenter = seq_fragmenter
        self.histogram = histogram
        self.rand_seed = rand_seed
        self.iter_count = 0

    def __iter__(self):
        self.iter_count += 1
        rng = np.random.RandomState(self.rand_seed)
        for seq in self.seq_generator.generator(rng):
            seq = self.mapper.apply(rng, seq)
            acgt_seq_splits = list(self.seq_fragmenter.get_acgt_seqs(seq))
            self.logger.debug('Splits of len={} to: {}'.format(len(seq), [len(f) for f in acgt_seq_splits]))

            for acgt_seq in acgt_seq_splits:
                kmer_seq = self.kmer_fragmenter.apply(rng, acgt_seq)  # list of strings
                if len(kmer_seq) > 0:
                    if self.iter_count == 1:
                        # only collect stats on the first call
                        self.histogram.add(kmer_seq)
                    yield kmer_seq


# In[51]:


### Histogram ###
import numpy as np
import os
from Bio.Seq import Seq

NUCLEOTIDES = 'ACGT'

class Tuple4:
    def __init__(self, pos1, pos2, neg1, neg2):
        self.pos1 = pos1
        self.pos2 = pos2
        self.neg1 = neg1
        self.neg2 = neg2

def determine_out_filename(output_dir, fileroot, mode, extension='txt'):
    return os.path.join(output_dir, '{}.{}.{}'.format(fileroot, mode, extension))

def create_tuple4(kmer1, kmer2, kmer1_neg, kmer2_neg):
    """
    all inputs are list of single nucleotides, e.g. ['A', 'A', 'C']
    """
    return Tuple4(
        ''.join(kmer1),
        ''.join(kmer2),
        ''.join(kmer1_neg),
        ''.join(kmer2_neg))

def insert_snippet(seq, snippet, idx):
    """
    idx: 0 <= idx <= len(seq)]
    """
    split1 = seq[:idx]
    split2 = seq[idx:]
    return split1 + snippet + split2

def pairwise_key(v1, v2):
    return '{}:{}'.format(v1, v2)

def rand_kmer(rng, k_low, k_high=None):
    """
    k_low and k_high are inclusive
    """
    if k_high is None:
        k_high = k_low
    k_len = rng.randint(k_low, k_high + 1)
    return ''.join([NUCLEOTIDES[x] for x in rng.randint(4, size=k_len)])

def rand_nt(rng):
    return NUCLEOTIDES[rng.randint(4)]

def generate_revcompl_pair(k_low, k_high=None, rng=None):
    # TODO make params k_high and rng be required
    if k_high is None:
        k_high = k_low
    if rng is None:
        rng = np.random
    kmer = rand_kmer(rng, k_low, k_high)
    return (kmer, revcompl(kmer))

def revcompl(kmer):
    return str(Seq(kmer).reverse_complement())

def generate_1nt_mutation_4tuple(rng, k_len):
    kmer1 = list(rand_kmer(rng, k_len, k_len))
    kmer2 = list(rand_kmer(rng, k_len, k_len))

    idx = rng.randint(len(kmer1))
    original_nt = kmer1[idx]
    mutate_nt = rand_nt(rng)

    kmer1_neg = list(kmer1)
    kmer1_neg[idx] = mutate_nt

    kmer2[idx] = mutate_nt
    kmer2_neg = list(kmer2)
    kmer2_neg[idx] = original_nt

    return create_tuple4(kmer1, kmer2, kmer1_neg, kmer2_neg)


# In[61]:


### Multi_K Model ###
from __future__ import print_function

import logbook
import tempfile
import numpy as np

from gensim.models import word2vec
from gensim import matutils

class SingleKModel:
    def __init__(self, model):
        self.model = model
        self.vocab_lst = sorted(model.vocab.keys())

class MultiKModel:
    def __init__(self, filepath):
        self.aggregate = word2vec.Word2Vec.load_word2vec_format(filepath, binary=False)
        self.logger = logbook.Logger(self.__class__.__name__)
        vocab_lens = [len(vocab) for vocab in self.aggregate.vocab.keys()]
        voc_vec = word2vec.Word2Vec(vocab, min_count=1)
        self.k_low = min(vocab_lens)
        self.k_high = max(vocab_lens)
        self.vec_dim = self.aggregate.vector_size

        self.data = {}
        for k in range(self.k_low, self.k_high + 1):
            self.data[k] = self.separate_out_model(k)

    def model(self, k_len):
        """
        Use vector('ACGTA') when possible
        """
        return self.data[k_len].model

    def vector(self, vocab):
        return self.data[len(vocab)].model[vocab]

    def unitvec(self, vec):
        return matutils.unitvec(vec)

    def cosine_distance(self, vocab1, vocab2):
        return np.dot(self.unitvec(self.vector(vocab1)), self.unitvec(self.vector(vocab2)))

    def l2_norm(self, vocab):
        return np.linalg.norm(self.vector(vocab))

    def separate_out_model(self, k_len):
        vocabs = [vocab for vocab in self.aggregate.vocab.keys() if len(vocab) == k_len]
        if len(vocabs) != 4 ** k_len:
            self.logger.warn('Missing {}-mers: {} / {}'.format(k_len, len(vocabs), 4 ** k_len))

        header_str = '{} {}'.format(len(vocabs), self.vec_dim)
        with tempfile.NamedTemporaryFile(mode='w') as fptr:
            print(header_str, file=fptr)
            for vocab in vocabs:
                vec_str = ' '.join("%f" % val for val in self.aggregate[vocab])
                print('{} {}'.format(vocab, vec_str), file=fptr)
            fptr.flush()
            return SingleKModel(word2vec.Word2Vec.load_word2vec_format(fptr.name, binary=False))


# In[ ]:


### Modified Execution Script ###
import glob
import logbook
from logbook.compat import redirect_logging
import configargparse
import numpy as np
from Bio import SeqIO

class arguements:
    def __init__(self, kmer_fragmenter, vec_dim, rseed, rseed_trainset, inputs,
                 k_low, k_high, context, epochs, gensim_iters, out_dir, debug):
        self.kmer_fragmenter = kmer_fragmenter
        self.vec_dim = vec_dim
        self.rseed = rseed
        self.rseed_trainset = rseed_trainset
        self.inputs = inputs
        self.k_low = k_low
        self.k_high = k_high
        self.context = context
        self.epochs = epochs
        self.gensim_iters = gensim_iters
        self.out_dir = out_dir
        self.debug = debug
#Setting Arguements for code
'''
kmer_fragmenter = 'sliding', vec_dim = 12, rseed=7, rseed_trainset=123, inputs,
                 k_low=5, k_high=5, context=4, epochs=1, gensim_iters=1, out_dir='./w2vEncoding', debug=False
'''
args = arguements(kmer_fragmenter = 'sliding', inputs = '../data/hg38/chroms/chr*.fa',k_low = 1,
                  k_high = 10, vec_dim = 100, epochs = 15, context = 10,
                  out_dir='../data/hg38/chroms', rseed = 7, rseed_trainset=123,
                 gensim_iters = 1, debug = False)

from gensim.models import word2vec

class Learner:
    def __init__(self, out_fileroot, context_halfsize, gensim_iters, vec_dim):
        self.logger = logbook.Logger(self.__class__.__name__)
        assert(word2vec.FAST_VERSION >= 0)
        self.logger.info('word2vec.FAST_VERSION (should be >= 0): {}'.format(word2vec.FAST_VERSION))
        self.model = None
        self.out_fileroot = out_fileroot
        self.context_halfsize = context_halfsize
        self.gensim_iters = gensim_iters
        self.use_cbow = 1
        self.vec_dim = vec_dim

        self.logger.info('Context window half size: {}'.format(self.context_halfsize))
        self.logger.info('Use CBOW: {}'.format(self.use_cbow))
        self.logger.info('gensim_iters: {}'.format(self.gensim_iters))
        self.logger.info('vec_dim: {}'.format(self.vec_dim))

    def train(self, kmer_seq_generator):
        self.model = word2vec.Word2Vec(
            sentences=kmer_seq_generator,
            size=self.vec_dim,
            window=self.context_halfsize,
            min_count=5,
            workers=4,
            sg=self.use_cbow,
            iter=self.gensim_iters)

        # self.logger.info(model.vocab)

    def write_vec(self):
        out_filename = '{}.w2v'.format(self.out_fileroot)
        self.model.save_word2vec_format(out_filename, binary=False)

def run_main(args, inputs, out_fileroot):
    logbook.info(' '.join(sys.argv))
    if not args.debug:
        import logging
        logging.getLogger('gensim.models.word2vec').setLevel(logging.INFO)

    np.random.seed(args.rseed)

    benchmark = Benchmark()

    if args.kmer_fragmenter == 'disjoint':
        kmer_fragmenter = DisjointKmerFragmenter(args.k_low, args.k_high)
    elif args.kmer_fragmenter == 'sliding':
        kmer_fragmenter = SlidingKmerFragmenter(args.k_low, args.k_high)
    else:
        raise InvalidArgException('Invalid kmer fragmenter: {}'.format(args.kmer_fragmenter))

    logbook.info('kmer fragmenter: {}'.format(args.kmer_fragmenter))

    histogram = Histogram()
    kmer_seq_iterable = KmerSeqIterable(
        args.rseed_trainset,
        SeqGenerator(inputs, args.epochs),
        SeqMapper(),
        SeqFragmenter(),
        kmer_fragmenter,
        histogram,
    )

    learner = Learner(out_fileroot, args.context, args.gensim_iters, args.vec_dim)
    learner.train(kmer_seq_iterable)
    learner.write_vec()

    histogram.print_stat(sys.stdout)

    benchmark.print_time()

def main():
    if args.debug:
        out_dir = '/tmp'
        log_level = 'DEBUG'
    else:
        out_dir = args.out_dir
        log_level = 'INFO'

    inputs = glob.glob(args.inputs)

    mbytes = estimate_bytes(inputs) // (10 ** 6)
    out_fileroot = get_output_fileroot(
        out_dir,
        'dna2vec',
        'k{}to{}-{}d-{}c-{}Mbp-{}'.format(
            args.k_low,
            args.k_high,
            args.vec_dim,
            args.context,
            mbytes * args.epochs,  # total Mb including epochs
            args.kmer_fragmenter))

    out_txt_filename = '{}.txt'.format(out_fileroot)
    with open(out_txt_filename, 'w') as summary_fptr:
        with Tee(summary_fptr):
            logbook.StreamHandler(sys.stdout, level=log_level).push_application()
            redirect_logging()
            run_main(args, inputs, out_fileroot)

if __name__ == '__main__':
    main()


# In[ ]:




