__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - Life Language Processing"
__website__ = "https://llp.berkeley.edu/"

import sys
sys.path.append('../')

import pickle as pickle
import codecs
import fnmatch
import os
from multiprocessing import Pool
import numpy as np
import tqdm
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from biom import load_table
from scipy import sparse
#from bootstrapping.bootstrapping import BootStrapping

class FileUtility(object):
    def __init__(self):
        print('File utility object created..')

    @staticmethod
    def generate_metagenomic_stat(file_directory, file_extenstion, only_files=[], proc=4, mapping_label=None):
        '''
        TODO: multicore needs to be fixed
        :param file_directory:
        :param file_extenstion:
        :param only_files:
        :param proc:
        :param mapping_label: to generate stat for each label separately
        :return:
        '''
        fasta_files, filename_mapping = FileUtility.read_fasta_directory(file_directory, file_extenstion,only_files=only_files)
        print(str(len(fasta_files)), 'fasta files found in', file_directory)
        print('Generating the stats ')
        file_sizes = dict()

        for fasta_file in tqdm.tqdm(fasta_files):
            name, seqs_size=FileUtility.read_sequence_file_length(fasta_file)
            if not mapping_label:
                if 'fatsa files' not in file_sizes:
                    file_sizes['fatsa files'] = []
                mapping_label['fatsa files'].append(seqs_size)
            else:
                for k in mapping_label.keys():
                    if k in name:
                        if not mapping_label[k] in file_sizes:
                            file_sizes[mapping_label[k]] = []
                        file_sizes[mapping_label[k]].append(seqs_size)

        output_list = list()
        output_list.append('BIOCPE AUTO GENERATED STATISTICS')
        for label in file_sizes.keys():
            output_list.append(' '.join(['number of ', label, ':', str(len(file_sizes[label]))]))
            output_list.append(
                ' '.join(['average number of sequences in ', label, ':', str(np.mean(file_sizes[label]))]))
            output_list.append(
                ' '.join(['std of number of sequences in ', label, ':', str(np.std(file_sizes[label]))]))
        FileUtility.save_list(file_directory + '/metastat.txt', output_list)

    @staticmethod
    def create_fasta_file(file_address, corpus, label):
        seq_id_pairs = [('.'.join([str(idx + 1), label[idx]]), x) for idx, x in enumerate(corpus)]
        seq_recs = [SeqRecord(Seq(seq, generic_dna), id=id, description='') for id, seq in seq_id_pairs]
        SeqIO.write(seq_recs, file_address, "fasta")


    @staticmethod
    def read_sequence_file(file_name_sample):
        '''
        :param file_name_sample:
        :return:
        '''
        corpus = []
        if file_name_sample[-1] == 'q':
            for cur_record in SeqIO.parse(file_name_sample, "fastq"):
                corpus.append(str(cur_record.seq).lower())
        else:
            for cur_record in SeqIO.parse(file_name_sample, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        return file_name_sample.split('/')[-1], corpus

    @staticmethod
    def read_sequence_file_length(file_name_sample):
        '''
        :param file_name_sample:
        :return:
        '''
        corpus = []
        if file_name_sample[-1] == 'q':
            for cur_record in SeqIO.parse(file_name_sample, "fastq"):
                corpus.append(str(cur_record.seq).lower())
        else:
            for cur_record in SeqIO.parse(file_name_sample, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        return file_name_sample.split('/')[-1], len(corpus)

    def train_cpe_merging(self, vocab_size, sample_size, directory):

        pool = Pool(processes=self.num_p)
        fasta_sample_files = [[x, sample_size] for x in self.fasta_files]
        corpus = []
        for ky, v in tqdm.tqdm(pool.imap_unordered(self.get_corpus, fasta_sample_files, chunksize=5),
                               total=len(self.fasta_files)):
            corpus = corpus + v
        print('Corpus size for training CPE is ', len(corpus))
        train_cpe(corpus, directory + 'cpe_16s_' + str(vocab_size), vocab_size,
                  directory + 'cpe_16s_' + str(vocab_size) + '_freq')

    @staticmethod
    def read_fasta_directory(file_directory, file_extenstion, only_files=[]):
        '''
        :param file_directory:
        :param file_extenstion:
        :param only_files:
        :return: list of fasta files, and a dic to map file to index
        '''
        if len(only_files) > 0:
            fasta_files = [x for x in FileUtility.recursive_glob(file_directory, '*.' + file_extenstion) if
                           x.split('/')[-1] in only_files]
        else:
            fasta_files = [x for x in FileUtility.recursive_glob(file_directory, '*.' + file_extenstion)]

        fasta_files.sort()
        mapping = {v: k for k, v in enumerate(fasta_files)}
        return fasta_files, mapping

    @staticmethod
    def read_OTU_format(biom_file):
        '''
        return OTU content
        :param biom_file:
        :return:
        '''
        table = load_table(biom_file)
        X_otu = table.matrix_data
        OTU_ID_Mapping = {x.split('.')[1]: idx for idx, x in enumerate(list(table.ids()))}
        return X_otu, OTU_ID_Mapping

    @staticmethod
    def save_obj(filename, value):
        with open(filename + '.pickle', 'wb') as f:
            pickle.dump(value, f)

    @staticmethod
    def load_obj(filename):
        return pickle.load(open(filename, "rb"))

    @staticmethod
    def save_list(filename, list_names):
        f = codecs.open(filename, 'w', 'utf-8')
        for x in list_names:
            f.write(x + '\n')
        f.close()

    @staticmethod
    def load_list(filename):
        return [line.strip() for line in codecs.open(filename, 'r', 'utf-8').readlines()]

    @staticmethod
    def save_sparse_csr(filename, array):
        np.savez(filename, data=array.data, indices=array.indices,
                 indptr=array.indptr, shape=array.shape)

    @staticmethod
    def load_sparse_csr(filename):
        loader = np.load(filename)
        return sparse.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])

    @staticmethod
    def _float_or_zero(value):
        try:
            return float(value)
        except:
            return 0.0

    @staticmethod
    def recursive_glob(treeroot, pattern):
        '''
        :param treeroot: the path to the directory
        :param pattern:  the pattern of files
        :return:
        '''
        results = []
        for base, dirs, files in os.walk(treeroot):
            good_files = fnmatch.filter(files, pattern)
            results.extend(os.path.join(base, f) for f in good_files)
        return results

    @staticmethod
    def read_fasta_sequences(file_name):
        corpus=[]
        for cur_record in SeqIO.parse(file_name, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        return corpus

    @staticmethod
    def read_fasta_sequences_ids(file_name):
        corpus=dict()
        for cur_record in SeqIO.parse(file_name, "fasta"):
                corpus[str(cur_record.id)]=(str(cur_record.seq).lower(),str(cur_record.description))
        return corpus
