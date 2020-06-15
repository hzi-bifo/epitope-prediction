from Bio import SeqIO
from pydpi.pypro import PyPro
from make_representations.sequencelist_representation import SequenceKmerRep, SequenceKmerEmbRep
from sklearn import svm, preprocessing
import sys
import numpy as np

protein = PyPro()


def readAAP(file):  #read AAP features from the AAP textfile
    aapdic = {}
    aapdata = open(file, 'r')
    for l in aapdata.readlines():
        aapdic[l.split()[0]] = float(l.split()[1])
    aapdata.close()
    return aapdic


def readAAT(file):  #read AAT features from the AAT textfile
    aatdic = {}
    aatdata = open(file, 'r')
    for l in aatdata.readlines():
        aatdic[l.split()[0][0:3]] = float(l.split()[1])
    aatdata.close()
    return aatdic


def aap(pep, aapdic, avg):  #return AAP features for the peptides
    feature=[]
    for a in pep:
        print(a)
        if int(avg) == 0:
            score = []
            count = 0
            for i in range(0, len(a) - 1):
                try:
                    score.append(round(float(aapdic[a[i:i + 2]]), 4))
                    # score += float(aapdic[a[i:i + 3]])
                    count += 1
                except KeyError:
                    # print(a[i:i + 3])
                    score.append(float(-1))
                    # score += -1
                    count += 1
                    continue
            # averagescore = score / count
            feature.append(score)
        if int(avg) == 1:
            score = 0
            count = 0
            for i in range(0, len(a) - 1):
                try:
                    score += float(aapdic[a[i:i + 2]])
                    count += 1
                except KeyError:
                    score += -1
                    count += 1
                    continue
            if count != 0:
                averagescore = score / count
            else:
                averagescore = 0
            feature.append(round(float(averagescore), 4))
    return feature


def aat(pep, aatdic, avg):  #return AAT features for the peptides
    feature = []
    for a in pep:
        if int(avg) == 0:
            # print(a)
            score = []
            count = 0
            for i in range(0, len(a) - 2):
                try:
                    score.append(round(float(aatdic[a[i:i + 3]]), 4))
                    # score += float(aapdic[a[i:i + 3]])
                    count += 1
                except KeyError:
                    # print(a[i:i + 3])
                    score.append(float(-1))
                    # score += -1
                    count += 1
                    continue
            # averagescore = score / count
            feature.append(score)
        if int(avg) == 1:
            score = 0
            count = 0
            for i in range(0, len(a) - 2):
                try:
                    score += float(aatdic[a[i:i + 3]])
                    count += 1
                except KeyError:
                    score += -1
                    count += 1
                    continue
            # print(a, score)
            if count != 0:
                averagescore = score / count
            else:
                averagescore = 0
            feature.append(round(float(averagescore), 4))
    return feature


def CTD(pep):  #Chain-Transition-Ditribution feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        ctd = protein.GetCTD()
        feature.append(list(ctd.values()))
    return feature


def AAC(pep):
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        aac = protein.GetAAComp()
        feature.append(aac)
    return feature


def DPC(pep):
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        dpc = protein.GetDPComp()
        feature.append(list(dpc.values()))
    return feature


def kmer(pep, k):
    feature = SequenceKmerRep(pep, 'protein', k)
    return feature


def protvec(pep, k, file):
    feature = SequenceKmerEmbRep(file, pep, 'protein', k)
    return feature


def readseq(file):  #read the sequence from the fasta file
    sequence = SeqIO.read(file, "fasta")
    return(str(sequence.seq))


def peptides(seq):  #return peptides of length 20 from the sequence
    pep = []
    i=0
    while i <= len(seq):
        if i+20 > len(seq):
            pep.append(seq[i:len(seq)])
        else:
            pep.append(seq[i:i+20])
        i = i + 20
    print(pep)
    return pep


def combinefeature(pep):
    aapdic = readAAP("aap-general.txt.normal")
    aatdic = readAAT("aat-general.txt.normal")
    f_aap = np.array(aap(pep, aapdic, 1))
    print(f_aap)
    f_aat = np.array(aat(pep, aatdic, 1))
    print(f_aat)
    f_aac = np.array(AAC(pep))
    print(f_aac)
    f_kmer = np.array(kmer(pep, 4).toarray())
    print(f_kmer)
    f_protvec = np.array(protvec(pep, 4, './protvec/uniref_3M.vec').toarray())
    print(f_protvec)
    return np.column_stack(f_aat,f_aac,f_kmer,f_protvec)


def scoremodel(file, mlfile):
    sequence = readseq(file)
    pep = peptides(sequence)
    features = combinefeature(pep)
    model = pickle.load(open(mlfile, 'rb'))
    pred_prob = model.fit(features)
    return pep, pred_prob


if __name__ == "__main__":
    peptide_list, pred_probability = scoremodel("test.fasta", "svm-model.pickle" )
    print("List of predicted epitopes:")
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            print(peptide_list[i], pred_probability[i][1])



    






