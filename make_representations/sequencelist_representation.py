__author__ = "Ehsaneddin Asgari"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu or ehsaneddin.asgari@helmholtz-hzi.de"
__project__ = "LLP - Life Language Processing"
__website__ = "https://llp.berkeley.edu/"

from alignment.pairwise_align import PWAlignment
from sklearn.feature_extraction.text import TfidfVectorizer
import itertools
from gensim.models import KeyedVectors
from scipy.sparse import csc_matrix
import numpy as np
from utility.file_utility import FileUtility


class SequenceKmerRep(object):
    def __init__(self, sequences, seq_type, k_mer, restricted_kmer=True, use_idf=False, norm=None, delete_empty_col=False):
        '''
            Class constructor
        '''
        self.seq_type = seq_type
        self.sequences = sequences
        self.k_mer = k_mer
        self.restricted_kmer = restricted_kmer
        self.delete_empty_col = delete_empty_col
        self.set_tfidf_vect(use_idf=use_idf, norm=norm)
        self.set_tfidf_representation()

    def set_tfidf_vect(self, use_idf=False, norm=None):
        '''
        :param use_idf: whether to use idf or not
        :param norm: whether to normalize or not
        :return:
        '''
        if self.restricted_kmer:
            if self.seq_type == 'protein':
                self.vocab = [''.join(xs) for xs in itertools.product(
                    'arndcqeghilkmfpstwyvbzxuo', repeat=self.k_mer)]
            if self.seq_type == 'dna':
                self.vocab = [''.join(xs) for xs in itertools.product(
                    'atcg', repeat=self.k_mer)]
            if self.seq_type == 'rna':
                self.vocab = [''.join(xs) for xs in itertools.product(
                    'aucg', repeat=self.k_mer)]
            self.vocab.sort()
            self.vectorizer = TfidfVectorizer(use_idf=use_idf, vocabulary=self.vocab, analyzer='char',
                                              ngram_range=(
                                                  self.k_mer, self.k_mer),
                                              norm=norm, stop_words=[], lowercase=True, binary=False)
        else:
            self.vectorizer = TfidfVectorizer(use_idf=use_idf, analyzer='char', ngram_range=(self.k_mer, self.k_mer),
                                              norm=norm, stop_words=[], lowercase=True, binary=False)

    def set_tfidf_representation(self):
        '''
        set the representation
        :return:
        '''
        self.X = self.vectorizer.fit_transform(self.sequences)
        self.vocab = self.vectorizer.get_feature_names()
        if self.delete_empty_col:
            del_col = np.where(np.sum(self.X.toarray(), axis=0) == 0)[0]
            self.X = self.X.toarray()
            self.X = np.delete(self.X, del_col, axis=1)
            self.X = csc_matrix(self.X)
            self.vocab = [v for i, v in enumerate(
                self.vocab) if i not in list(del_col)]

    def get_representation(self, seq):
        '''
        :param seq: representation for an extra sequence
        :return:
        '''
        return self.vectorizer.fit_transform([seq]).toarray()


class SequenceKmerEmbRep(SequenceKmerRep):
    def __init__(self, embedding_file, sequences, seq_type, k_mer, restricted_kmer=True, use_idf=False, norm=None):
        '''
            Class constructor
        '''
        SequenceKmerRep.__init__(self, sequences, seq_type, k_mer, restricted_kmer=restricted_kmer, use_idf=use_idf,
                                 norm=norm, delete_empty_col=True)
        self.model = KeyedVectors.load_word2vec_format(
            embedding_file, binary=False)
        model = self.model
        try:
            k_mer_dict = FileUtility.load_obj(
                '../config/' + str(k_mer) + "_in_model")
        except:
            k_mer_dict = dict()
        new_words = [x.lower() for x in (self.vocab) if x.upper() not in model]
        for w in new_words:
            if w not in k_mer_dict:
                k_mer_dict[w] = self.closest_kmer_in_model(w)
        FileUtility.save_obj('../config/' + str(k_mer) +
                             "_in_model", k_mer_dict)
        # produce embedding mapping
        self.emb_trans = [self.model[x.upper()] if x.upper(
        ) in self.model else self.model[k_mer_dict[x]] for x in self.vocab]
        # summation vector
        self.embeddingX = self.X.dot(self.emb_trans)

    def closest_kmer_in_model(self, k_mer):
        '''
            Look for the closest k-mer
        :param k_mer:
        :return:
        '''
        k_mer = k_mer.upper()
        value = -1000
        closest = ''
        for w in self.model.index2word:
            if '<' not in w:
                PWA = PWAlignment(k_mer, w, self.seq_type)
                if PWA.getAlignScore() > value:
                    value = PWA.getAlignScore()
                    closest = w
        return closest
