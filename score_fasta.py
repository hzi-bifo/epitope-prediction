from Bio import SeqIO
from pydpi.pypro import PyPro
from make_representations.sequencelist_representation import SequenceKmerRep, SequenceKmerEmbRep
from sklearn.metrics import precision_score, recall_score, roc_auc_score, auc, matthews_corrcoef, classification_report
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, roc_curve, precision_recall_curve, precision_recall_fscore_support
from sklearn import svm, preprocessing
from sklearn.model_selection import KFold, GridSearchCV, StratifiedKFold
import sys
import numpy as np
import os.path
import pickle
import pylab as pl
from scipy import interp

protein = PyPro()

class MyCustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "__main__":
            module = "score_fasta"
        return super().find_class(module, name)




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
        feature.append(list(aac.values()))
    return feature


def DPC(pep): # Dipeptide Composition feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        dpc = protein.GetDPComp()
        feature.append(list(dpc.values()))
    return feature


def kmer(pep, k, testing, vocab): # Calculate k-mer feature
    feature = SequenceKmerRep(pep, 'protein', k, vocab=vocab, testing=testing)
    return feature


def protvec(pep, k, file): #Calculate ProtVec representation
    feature = SequenceKmerEmbRep(file, pep, 'protein', k)
    return feature
    
def PAAC(pep):
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        paac=protein.GetMoranAuto()
        #paac = protein.GetPAAC(lamda=4)
        feature.append(list(paac.values()))
        name = list(paac.keys())
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
    while i < len(seq):
        if len(seq)<20:
            pep.append(seq)
            break
        if len(seq)>=20:
            if i+20 > len(seq):
                break
            else:
                pep.append(seq[i:i+20])
            i = i + 1
    print(pep)
    return pep
    

def precision_0(y_true, y_pred, labels=None, average='binary', sample_weight=None):
    '''
    :param y_true:
    :param y_pred:
    :param labels:
    :param average:
    :param sample_weight:
    :return: calculate prec for neg class
    '''
    p, _, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 beta=1,
                                                 labels=labels,
                                                 pos_label=0,
                                                 average=average,
                                                 warn_for=('f-score',),
                                                 sample_weight=sample_weight)
    return p


def recall_0(y_true, y_pred, labels=None, average='binary', sample_weight=None):
    '''
    :param y_true:
    :param y_pred:
    :param labels:
    :param average:
    :param sample_weight:
    :return: calculate recall for neg class
    '''
    _, r, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 beta=1,
                                                 labels=labels,
                                                 pos_label=0,
                                                 average=average,
                                                 warn_for=('f-score',),
                                                 sample_weight=sample_weight)
    return r


def f1_0(y_true, y_pred, labels=None, average='binary', sample_weight=None):
    '''
    :param y_true:
    :param y_pred:
    :param labels:
    :param average:
    :param sample_weight:
    :return: calculate f1 for neg class
    '''
    _, _, f, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 beta=1,
                                                 labels=labels,
                                                 pos_label=0,
                                                 average=average,
                                                 warn_for=('f-score',),
                                                 sample_weight=sample_weight)
    return f




def readmodel(mlfile):
    '''try:
         print(mlfile)
         return pickle.load(open(mlfile, 'rb'))
    except:
        print("Error in reading model file")
        sys.exit()'''
    with open(mlfile, 'rb') as f:
        unpickler = MyCustomUnpickler(f)
        obj = unpickler.load()
    return obj
    

def combinefeature(pep, featurelist, vocab, aap_file, aat_file):
    print (featurelist)
    a=np.empty([len(pep), 1])
    if 'aap' in featurelist:
        aapdic = readAAP("./aap/"+aap_file)
        f_aap = np.array([aap(pep, aapdic, 1)]).T
        a = np.column_stack((a,f_aap))
        #print(f_aap)
    if 'aat' in featurelist:
        aatdic = readAAT("./aat/"+aat_file)
        f_aat = np.array([aat(pep, aatdic, 1)]).T
        a = np.column_stack((a, f_aat))
        #print(f_aat)
    if 'dpc' in featurelist:
        f_dpc, name = DPC(pep)
        # f_dpc = np.average(f_dpc, axis =1)
        a = np.column_stack((a, np.array(f_dpc)))
    if 'aac' in featurelist:
        f_aac = AAC(pep)
        a = np.column_stack((a, np.array(f_aac)))
        #fname = fname + name
    if 'paac' in featurelist:
        f_paac = PAAC(pep)
        f_paac = pca.fit_transform(f_paac)
        a = np.column_stack((a, np.array(f_paac)))
        #fname = fname + name
    if 'kmer' in featurelist:
        kmers = kmer(pep, 4, vocab=vocab, testing=1)
        #f_kmer = np.array(kmers.X.toarray())
        f_kmer = np.array(kmers.X.toarray())
        a = np.column_stack((a, f_kmer))
        #fname = fname + name
    if 'ctd' in featurelist:
        f_ctd, name = CTD(pep)
        a = np.column_stack((a, np.array(f_ctd)))
        #fname = fname + name
    if 'protvec' in featurelist:
            if os.path.isfile('./protvec/sp_sequences_4mers_vec.bin') == True:
                f_protvec = np.array(protvec(pep, 4, './protvec/sp_sequences_4mers_vec.bin').embeddingX)
                a = np.column_stack((a, f_protvec))
            else:
                print("Protvec binaries are missing. See README file.")
                sys.exit()    
                #print(f_protvec)
    #print(a)
    return a[:,1:]


def plot(model, x,y):
    cv = StratifiedKFold(n_splits=5)
    splits = list(cv.split(x,y))
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []
    classifier = model
    for i, (train, test) in enumerate(splits):
        probas_ = classifier.fit(x[train], y[train]).predict_proba(x[test])
        #clf.fit(x[train], y[train])
        predict_values = classifier.fit(x[train], y[train]).predict(x[test])
        print(classification_report(y[test], predict_values))
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        pl.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))
    pl.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
    mean_tpr /= len(splits)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    pl.plot(mean_fpr, mean_tpr, 'k--',
            label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

    pl.xlim([-0.05, 1.05])
    pl.ylim([-0.05, 1.05])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('BCPred')
    pl.show()
    

def predict(training_data, features):
    model = training_data['model'].best_estimator_
    training_features = training_data['training_features']
    scaling = training_data['scaling']
    #plot(model , scaling.transform(training_features), np.array([1]*701 + [0]*701))
    features = scaling.transform(features)
    #print(model.score(scaling.transform(training_features),([1]*701)+([0]*701)))
    try:
        return model.predict_proba(features)
    except:
        print("Error in predicting epitopes.")
        sys.exit()
        
def test(training_data, x_test, y_test):
    model = training_data['model'].best_estimator_
    scaling = training_data['scaling']
    features = scaling.transform(x_test)
    print (model.score(features, y_test))


def scoremodel(file, mlfile, aap_file, aat_file):
    sequence = readseq(file)
    pep = peptides(sequence)
    training_data= readmodel(mlfile)
    #print(training_data.keys())
    features = combinefeature(pep, training_data['featurelist'], training_data['vocab'], aap_file, aat_file)
    '''newdata = open('abcpred-20.txt', 'r')
    anew = []
    for l in newdata.readlines():
        if l[0] == '#':
            continue
        else:
            anew.append(l.strip().split())
    anew = np.array(anew)
    newdata.close()
    testpeptides = anew[:,0]
    y_test = anew[:, -1].astype(int)
    x = combinefeature(testpeptides, training_data['featurelist'], training_data['vocab'])
    print(len(x),len(y_test))
    test(training_data, x, y_test)'''
    return pep, predict(training_data, features)


if __name__ == "__main__":
    peptide_list, pred_probability = scoremodel("./input/example.fasta", "./model/svm-model.pickle" )
    print("List of predicted epitopes:")
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            print(peptide_list[i], pred_probability[i][1])
