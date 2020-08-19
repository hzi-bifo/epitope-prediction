from Bio import SeqIO
from pydpi.pypro import PyPro
from make_representations.sequencelist_representation import SequenceKmerRep, SequenceKmerEmbRep
from sklearn import svm, preprocessing
import sys
import numpy as np
import os.path

protein = PyPro()


def readAAP(file):  #read AAP features from the AAP textfile
    try:
        aapdic = {}
        aapdata = open(file, 'r')
        for l in aapdata.readlines():
            aapdic[l.split()[0]] = float(l.split()[1])
        aapdata.close()
        return aapdic
    except:
        print("Error in reading AAP feature file. Please make sure that the AAP file is correctly formatted")
        sys.exit()


def readAAT(file):  #read AAT features from the AAT textfile
    try:
        aatdic = {}
        aatdata = open(file, 'r')
        for l in aatdata.readlines():
            aatdic[l.split()[0][0:3]] = float(l.split()[1])
        aatdata.close()
        return aatdic
    except:
        print("Error in reading AAT feature file. Please make sure that the AAT file is correctly formatted")
        sys.exit()


def aap(pep, aapdic, avg):  #return AAP features for the peptides
    feature=[]
    for a in pep:
        #print(a)
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


def AAC(pep): # Single Amino Acid Composition feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        aac = protein.GetAAComp()
        feature.append(aac)
    return feature


def DPC(pep): # Dipeptide Composition feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        dpc = protein.GetDPComp()
        feature.append(list(dpc.values()))
    return feature


def kmer(pep, k): # Calculate k-mer feature
    feature = SequenceKmerRep(pep, 'protein', k)
    return feature


def protvec(pep, k, file): #Calculate ProtVec representation
    feature = SequenceKmerEmbRep(file, pep, 'protein', k)
    return feature


def readseq(file):  #read the sequence from the fasta file
    try:
        sequence = SeqIO.read(file, "fasta")
        for i in sequence.seq:
            #print(i)
            if i in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] :
                continue
            else:
                print("Invalid amino acid code found. Please enter sequences with only 20 aa code.")
                sys.exit()
        return(str(sequence.seq))
    except ValueError:
        print("Please enter a valid fasta file")
        sys.exit()


def peptides(seq):  #return peptides of length 20 from the sequence
    pep = []
    i=0
    while i <= len(seq):
        if i+20 > len(seq):
            pep.append(seq[i:len(seq)])
        else:
            pep.append(seq[i:i+20])
        i = i + 20
    #print(pep)
    return pep


def readmodel(mlfile):
    try:
        return pickle.load(open(mlfile, 'rb'))
    except:
        print("Error in reading model file")
        sys.exit()


def predict(model, features):
    return model.predict(features)
    '''except:
        print("Error in predicting epitopes.")
        sys.exit()'''


def combinefeature(pep):
    aapdic = readAAP("./aap/aap-general.txt.normal")
    aatdic = readAAT("./aat/aat-general.txt.normal")
    f_aap = np.array(aap(pep, aapdic, 1))
    #print(f_aap)
    f_aat = np.array(aat(pep, aatdic, 1))
    #print(f_aat)
    f_aac = np.array(AAC(pep))
    #print(f_aac)
    f_kmer = np.array(kmer(pep, 4).X.toarray())
    #print(f_kmer)
    if os.path.isfile('./protvec/sp_sequences_4mers_vec.txt') == True:
        f_protvec = np.array(protvec(pep, 4, './protvec/sp_sequences_4mers_vec.txt').embeddingX)
    else:
        print("Protvec binaries are missing. See README file.")
        sys.exit()    
        #print(f_protvec)
    return np.column_stack((f_aat,f_aac,f_kmer,f_protvec))


def scoremodel(file, mlfile):
    sequence = readseq(file)
    pep = peptides(sequence)
    features = combinefeature(pep)
    #print(len(features[0]))
    model = readmodel(mlfile)
    return pep, predict(model, features)


if __name__ == "__main__":
    peptide_list, pred_probability = scoremodel("./input/example.fasta", "./model/svm-model.pickle" )
    print("List of predicted epitopes:")
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            print(peptide_list[i], pred_probability[i][1])