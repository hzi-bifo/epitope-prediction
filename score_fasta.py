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
#import pylab as pl
from scipy import interp

protein = PyPro()

class MyCustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "__main__":
            module = "score_fasta"
        return super().find_class(module, name)


def error(err):
    sys.exit(err)


def readAAP(file):  #read AAP features from the AAP textfile
    try:
        aapdic = {}
        aapdata = open(file, 'r')
        for l in aapdata.readlines():
            aapdic[l.split()[0]] = float(l.split()[1])
        aapdata.close()
        return aapdic
    except:
        error("Error in reading AAP feature file. Please make sure that the AAP file is correctly formatted")


def readAAT(file):  #read AAT features from the AAT textfile
    try:
        aatdic = {}
        aatdata = open(file, 'r')
        for l in aatdata.readlines():
            aatdic[l.split()[0][0:3]] = float(l.split()[1])
        aatdata.close()
        return aatdic
    except:
        error("Error in reading AAT feature file. Please make sure that the AAT file is correctly formatted")


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


def combinefeature(pep, featurelist, vocab, aapdic, aatdic):
    
    a=np.empty([len(pep), 1])
    if 'aap' in featurelist:
        #aapdic = readAAP("./aap/"+aap_file)
        f_aap = np.array([aap(pep, aapdic, 1)]).T
        a = np.column_stack((a,f_aap))
        #print(f_aap)
    if 'aat' in featurelist:
        #aatdic = readAAT("./aat/"+aat_file)
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
        f_protvec = np.array(protvec(pep, 4, './protvec/sp_sequences_4mers_vec.bin').embeddingX)
        a = np.column_stack((a, f_protvec))
            
    #print(a)
    return a[:,1:]


def readseq(file, length):  #read the sequences from the fasta file
    print("=================================================================\n")
    try:
        records = list(SeqIO.parse(file, "fasta"))
        # print("Amino acid sequence:",sequence.seq)
        sequences={}
        
        if len(records) > 1:
            print("Found "+str(len(records))+" protein sequences\n")
        else:
            print("Found "+str(len(records))+" protein sequence\n")
        
        maxlenname=0
        maxlenseq=0
        for i in records:
            if len(str(i.id)) > maxlenname:
                maxlenname = len(str(i.id))
            if len(str(i.seq)) > maxlenseq:
                maxlenseq = len(str(i.seq))
        
        for i in records:
            #print("test")
            seq_check = 0
            len_check = 0
            if (len(str(i.seq)) >=4):
                len_check = 1
            if (len(str(i.seq)) <4):
                len_check = 0
                print(i.id.ljust(maxlenname,' '),":",i.seq)
                print("Error in "+i.id+" sequence.The length of the sequence can't be less than 4. This sequence will be skipped. \n", file=sys.stderr)
                if(len(records))==1:
                        error("Error. No valid sequence found")
                continue
            for j in i.seq:
                if j.upper() in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] :
                    seq_check = 1
                else:
                    seq_check = 0
                    print(i.id.ljust(maxlenname,' '),":",i.seq)
                    print("Error in "+i.id+" sequence. \""+j+"\" is not a valid amino acid. This sequence will be skipped. Please enter sequences with only 20 aa code. \n", file=sys.stderr)
                    if(len(records))==1:
                        error("Error. No valid sequence found")
                    break
            if seq_check == 1 and len_check ==1 :
                sequences[i.id]=i.seq
                print(i.id.ljust(maxlenname,' '),":",i.seq,"\n")
        print("Number of valid sequences found:",len(sequences.keys()),"\n")
        print("=================================================================\n")
        return(sequences, maxlenname, maxlenseq)
    except ValueError:
        error("Error in reading the input fasta file. Please enter a valid fasta file")




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
    try:
         print("Loading modelfile:", mlfile,"\n")
         with open(mlfile, 'rb') as f:
             unpickler = MyCustomUnpickler(f)
             obj = unpickler.load()
             return obj
    except:
        error("Error in reading model file")


def calcpeptides(seq, length, shortcheck, longcheck):  #return peptides of length 20 from the sequence
    pep = []
    peploc=[]
    i=0
    while i < len(seq):
        if len(seq)<length:
            shortcheck = 1
            pep.append(seq.upper())
            peploc.append(i+1)
            break
        if len(seq)>=length:
            if i+length > len(seq):
                break
            else:
                pep.append(seq[i:i+length].upper())
                peploc.append(i+1)
                longcheck = 1
            i = i + 1
    #print(pep)
    return pep, peploc, shortcheck, longcheck
    


def createpeptides(sequences, length, maxname, maxlen): #calculates petides for all the sequences
    pep = []
    peploc = []
    seqid = []
    print("Generating peptides of length "+str(length)+"\n")
    global_longcheck = 0
    global_shortcheck = 0
    digit_len=len(str(maxlen))
    digit_pepno=len(str((maxlen-length+1)))
    for i in sequences.keys():
        shortcheck = 0
        longcheck = 0
        peptemp, peploctemp, shortcheck, longcheck  = calcpeptides(str(sequences[i]), length, shortcheck, longcheck)
        pep = pep + peptemp
        peploc = peploc + peploctemp
        seqid = seqid + [i]*len(peptemp)
        seq_len = len(str(sequences[i]))
        no_pep = len(peptemp)
        if longcheck == 1 and shortcheck == 0:
            print("sequence:",i.center(maxname,' '),"| length: "+"{:{align}{width}}".format(str(seq_len),align='>', width = digit_len),"| no. of peptides: "+"{:{align}{width}}".format(str(no_pep),align='<', width = digit_pepno)+"\n")
            global_longcheck = 1
        if shortcheck == 1:
            print("sequence:",i.center(maxname,' '),"| length: "+"{:{align}{width}}".format(str(seq_len),align='>', width = digit_len),"| length is smaller than the specified peptide length of "+str(length)+"\n")
            global_shortcheck = 1



    print()
    if global_longcheck == 0:
        print("No sequences with length >",length,"found. Predictions will be made on the entire sequences\n")

    if global_shortcheck == 1:
        print("Note: If the length of a sequence is smaller than the specified length of the peptide (in this case: "+str(length)+"), the prediction will be made on the entire sequence.\n")
    print("=================================================================\n")
    return pep, peploc, seqid
    
    
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
    return model.predict_proba(features)
    #print(model.score(scaling.transform(training_features),([1]*701)+([0]*701)))
    '''try:
        return model.predict_proba(features)
    except:
        print("Error in predicting epitopes.")
        sys.exit()'''
        

def test(training_data, x_test, y_test):
    model = training_data['model'].best_estimator_
    scaling = training_data['scaling']
    features = scaling.transform(x_test)
    print (model.score(features, y_test))


def scoremodel(file, mlfile, length):
    sequences, maxname, maxlen = readseq(file, length)
    pep, peploc, seqid = createpeptides(sequences, length, maxname, maxlen)
    training_data= readmodel(mlfile)
    #print(training_data.keys())
    print ("Calculating the",training_data['featurelist'],"features for the peptides. \n")
    #print(*training_data['featurelist'], sep=",")
    features = combinefeature(pep, training_data['featurelist'], training_data['vocab'],training_data['aap'],training_data['aat'])
    print("\nFeature calculation complete!\n")
    print("=================================================================\n")
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
    print("Making predictions \n")
    return sequences, pep, peploc, seqid, predict(training_data, features)


def print_stdout(pep, peploc, seqid, pred, cutoff, length):
    epitopelist = []
    nonepitopelist = []
    for i in range(len(pred)):
        if pred[i][1] >= cutoff:
            epitopelist.append(i)
        if pred[i][1] < cutoff:
            nonepitopelist.append(i)
    
    if len(epitopelist) > 0:
        print("Peptides predicted as epitopes: \n")
        #print("seq_id"+"\t"+"peptide"+"\t"+"start"+"\t"+"end"+"\t"+"score")
        seqid_maxname = 0
        pep_maxlen = 0
        startloc_max = 0
        for i in epitopelist:
            if len(str(seqid[i])) > seqid_maxname:
                seqid_maxname = len(str(seqid[i]))
            if len(str(pep[i])) > pep_maxlen:
                pep_maxlen = len(str(pep[i]))
            if peploc[i] > startloc_max:
                startloc_max = peploc[i]
        width_seqid = 6
        if seqid_maxname > 6:
            width_seqid = seqid_maxname
        width_pep = pep_maxlen
        width_startloc = 5
        if len(str(startloc_max)) > 5:
            width_startloc = len(str(startloc_max))
        width_endloc = 5
        if len(str(startloc_max+length-1)) > 5:
            width_endloc = len(str(startloc_max+length-1))
        print("seq_id".rjust(width_seqid,' ')+"\t"+"peptide".center(width_pep,' ')+"\t"+"start".rjust(width_startloc,' ')+"\t"+"end".rjust(width_endloc,' ')+"\t"+"score".ljust(5,' '))
        print("".rjust(width_seqid,'-')+"\t"+"".center(width_pep,'-')+"\t"+"".rjust(width_startloc,'-')+"\t"+"".rjust(width_endloc,'-')+"\t"+"".ljust(5,'-'))
        for i in epitopelist:
            print(str(seqid[i]).rjust(width_seqid,' ')+"\t"+str(pep[i]).center(width_pep,' ')+"\t"+str(peploc[i]).rjust(width_startloc,' ')+"\t"+str(peploc[i]+len(str(pep[i]))-1).rjust(width_endloc,' ')+"\t"+'{:{width}.{prec}f}'.format(pred[i][1], width = 5, prec = 3))
        print()
    else:
        print("No peptides were predicted as epitopes \n")


def print_file(sequences, pep, peploc, seqid, pred, cutoff, length, output):
    #print("=================================================================\n")
    #print("Writing to files.")
    epitopelist = []
    for i in range(len(pred)):
        if pred[i][1] >= cutoff:
            epitopelist.append(i)

    seqid_maxname=0
    for i in sequences.keys():
            if len(str(i)) > seqid_maxname:
                seqid_maxname = len(str(i))




    if len(epitopelist) > 0:
        output_string={}
        for i in sequences.keys():
            output_string[i]=['',0]
        for i in epitopelist:
            output_string[seqid[i]][0] = output_string[seqid[i]][0]+(str(seqid[i])+"\t"+str(pep[i])+"\t"+str(peploc[i])+"\t"+str(peploc[i]+len(str(pep[i]))-1)+"\t"+'{:{width}.{prec}f}'.format(pred[i][1], width = 5, prec = 3)+"\n")
            output_string[seqid[i]][1] = output_string[seqid[i]][1] + 1
    else:
        print("None of the peptides were predicted as epitopes. \n")
        sys.exit()

    for i in output_string.keys():
        with open(output+"/"+str(i)+".epitopes.txt", "w") as fileoutput:
            if output_string[i][0] != '':
                print("sequence:",i.center(seqid_maxname,' '),"| no. of peptides predicted as epitopes:",output_string[i][1])
                print("writing predictions to:",output+"/"+str(i)+".epitopes.txt\n")
                fileoutput.write("seq_id"+"\t"+"peptide"+"\t"+"start"+"\t"+"end"+"\t"+"score\n")
                fileoutput.write(output_string[i][0])
            else:
                print("sequence:",i.center(seqid_maxname,' '),"| no. of peptides predicted as epitopes:",output_string[i][1])
                print("writing predictions to:",output+"/"+str(i)+".epitopes.txt\n")
                fileoutput.write("None of the peptides were predicted as epitopes.")






