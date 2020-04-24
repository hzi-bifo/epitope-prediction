import warnings
warnings.filterwarnings("ignore")
from sklearn import svm, datasets, metrics
from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, roc_curve, precision_recall_curve, precision_recall_fscore_support, roc_auc_score
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import classification_report
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold, GridSearchCV, StratifiedKFold
# from sklearn.cross_validation import 
import numpy as np
import sys
import pylab as pl
import matplotlib.pyplot as plt
#from sklearn.grid_search import GridSearchCV
stdsc = StandardScaler()
from matplotlib.colors import ListedColormap
from sklearn.neural_network import MLPClassifier
from scipy import interp
from sklearn.metrics import matthews_corrcoef as mcc
from pydpi.pypro import PyPro
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
# from classifier.classical_classifiers import RFClassifier, SVM
# from make_representations.sequencelist_representation import SequenceKmerRep, SequenceKmerEmbRep
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import precision_score, recall_score
import argparse


protein = PyPro()

parser = argparse.ArgumentParser(
    description=("Take training data and give scores "
                 "can take additonal test data"),
    formatter_class=argparse.RawTextHelpFormatter,
    usage=("%(prog)s inputfile [-h]"
           "[-i training data][-aat AAT][-r REDUCED] [-a Average/extended] [-aap AAP] [-t testing data"))
parser.add_argument("-i", "--input", help=(" input training file"))
parser.add_argument("-t", "--test", help=(" test data"))
parser.add_argument("-aap", "--AAP", help=(" AAP scale"))
parser.add_argument("-aat", "--AAT", help=(" AAT scale"))
parser.add_argument("-r", "--reduce", help=(" reduce scales"))
parser.add_argument("-a", "--average", help=(" averaged or extended feature"))
args = parser.parse_args()

'''relevant=open('temp.txt','r')
rel=[]
for l in relevant.readlines():
    rel.append(l.strip())
print(rel)'''

aapdic = {}
aapdata = open(args.AAP, 'r')
for l in aapdata.readlines():
    aapdic[l.split()[0]] = float(l.split()[1])
aapdata.close()

aatdic = {}
aatdata = open(args.AAT, 'r')
for l in aatdata.readlines():
    aatdic[l.split()[0][0:3]] = float(l.split()[1])
aatdata.close()


def calculate_aat_scale(poslines, neglines, reduce):
    poscount = {}
    negcount = {}
    for i in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
        for j in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            for k in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
                poscount[i + j + k] = 1
                negcount[i + j + k] = 1
    postotal = 8000
    for l in poslines:
        seq = l.strip()
        for i in range(0, len(l) - 3):
            tp = seq[i: i + 3].upper()
            if 'X' in tp.upper():
                continue
            try:
                poscount[tp] = poscount[tp] + 1
                postotal = postotal + 1
            except KeyError:
                continue
    print(poscount)
    negtotal = 8000
    seq = ''
    seqcount = 0
    for line in range(1, len(neglines)):
        seq = seq + neglines[line].strip()
        if neglines[line][0] == '>' or (line == len(neglines) - 1):
            seqcount += 1
            seq = seq.split('>')[0].upper()
            for i in range(0, len(seq) - 3):
                tp = seq[i: i + 3]
                if 'X' in tp.upper():
                    continue
                try:
                    negcount[tp] = negcount[tp] + 1
                    negtotal = negtotal + 1
                except KeyError:
                    continue
            seq = ''
            continue
    aatscale = {}
    for i in poscount.keys():
        try:
            aatscale[i] = (poscount[i] / postotal) / (negcount[i] / negtotal)
        except KeyError:
            continue

    return aatscale


def aap(a, aapdic, avg, r):
    # print(aapdic)
    if int(r) == 1:
        group = {'R': 'B', 'K': 'B', 'E': 'J', 'D': 'J', 'S': 'O',
                 'T': 'O', 'L': 'U', 'V': 'U', 'I': 'U', 'P': 'P',
                 'Q': 'X', 'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A',
                 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'Y': 'Y',
                 'X': 'X'}
        seq = a
        a = ''
        for i in seq:
            try:
                a += group[i]
            except KeyError:
                a += i
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
        return score
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
        return round(float(averagescore), 4)


def aat(a, aatdic, avg, r):
    if int(r) == 1:
        group = {'R': 'B', 'K': 'B', 'E': 'J', 'D': 'J', 'S': 'O',
                 'T': 'O', 'L': 'U', 'V': 'U', 'I': 'U', 'P': 'P',
                 'Q': 'X', 'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A',
                 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'Y': 'Y',
                 'X': 'X'}
        seq = a
        a = ''
        for i in seq:
            try:
                a += group[i]
            except KeyError:
                a += i
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
        return score
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
        return round(float(averagescore), 4)


def writedata(a):
    a = np.array(a)
    x = np.array(preprocessing.scale(a[:, 1:-1]))
    outfile = open("training-set" + ".feat", 'w')
    for i in range(len(a)):
        print(a[i, 0] + "\t", end='', file=outfile)
        for j in range(len(x[0])):
            print(str(x[i, j]) + "\t", end='', file=outfile)
        print(a[i, -1] + "\t", end='', file=outfile)
        print(file=outfile)
    outfile.close()


def CTD(seq):
    protein.ReadProteinSequence(seq)
    ctd = protein.GetCTD()
    return list(ctd.values())


def AAC(seq):
    protein.ReadProteinSequence(seq)
    aac = protein.GetAAComp()
    return aac


def DPC(seq):
    protein.ReadProteinSequence(seq)
    dpc = protein.GetDPComp()
    return list(dpc.values())


data = open(args.input, 'r')
a = []
for l in data.readlines():
    if l[0] == '#':
        continue
    else:
        a.append(l.strip().split())

x = []
y = []
# print(a)
# print (a[0])
# make a function for making features
for i in range(len(a)):
    if (args.AAP != 0):
        aapscore = aap(a[i][0], aapdic, args.average, args.reduce)
        a[i] = np.insert(a[i], len(a[i]) - 1, aapscore)
    if args.AAT != 0:
        aatscore = aat(a[i][0], aatdic, args.average, args.reduce)
        a[i] = np.insert(a[i], len(a[i]) - 1, aatscore)
    # ctd = CTD(a[i][0])
    # a[i] = np.insert(a[i], len(a[i]) - 1, ctd)

    '''aac = [AAC(a[i][0])[x]
           for x in ['A', 'F', 'G', 'I', 'L', 'N', 'P', 'Q', 'T', 'W']]
    a[i] = np.insert(a[i], len(a[i]) - 1, aac)

    dpc = DPC(a[i][0])
    a[i] = np.insert(a[i], len(a[i]) - 1, dpc)'''

a = np.array(a)
print("test")
print(a)
y = a[:, -1].astype(int)
x = a[:, 1:-1].astype(float)
seq = a [:, 0]
# sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
# x = x.reshape(-1,1)

# print(x)
y = np.array(y)
writedata(a)

'''clf = ExtraTreesClassifier(n_estimators=50)
clf = clf.fit(x, y)
model = SelectFromModel(clf, prefit=True)
x = model.transform(x)
#x = SelectKBest(chi2, k=10).fit_transform(x, y)
# print(x[0])
# x = sel.fit_transform(x)
# print(x[0])'''
x = np.array(preprocessing.scale(x))

# save pickle file
# numpy.save(file, arr, allow_pickle=True, fix_imports=True)[source]


def svc_param_selection(x, y, cv):
    ''' Cs = np.logspace(-2, 10, 13)
    gammas = np.logspace(-9, 3, 13)
    param_grid = {'C': Cs, 'gamma': gammas}'''
    grid_search = GridSearchCV(SVC(kernel='rbf', probability=True),
                               param_grid={'C': [1000, 500, 200, 100, 50,
                                                 20, 10, 2, 1, 0.2, 0.5,
                                                 0.01, 0.02, 0.05, 0.001],
                                           'gamma': [1000, 500, 200, 100,
                                                     50, 20, 10, 5, 2, 1,
                                                     0.2, 0.5, 0.01, 0.02,
                                                     0.05, 0.001, 0.0001]},
                               scoring='accuracy', cv=cv, n_jobs=8)
    grid_search.fit(x, y)
    # grid_search.best_params_
    return grid_search.best_params_, grid_search.best_score_


cv = StratifiedKFold(n_splits=5)
splits = list(cv.split(x,y))
# print(svc_param_selection(x, y, cv))
classifier = svm.SVC(kernel='rbf', C=1000, gamma=0.001, probability=True)
# classifier = svm.SVC(kernel='linear', C=1000, gamma=0.0001, probability=True)
sequences = a[:, 0]
# SKR = SequenceKmerRep(sequences, 'protein', 4)

mean_tpr = 0.0
mean_fpr = np.linspace(0, 1, 100)
all_tpr = []

print("Number of Features:", len(x[0]))
print(x, y)

score = 0
mccscore = 0
scorenew = 0
mccscorenew = 0
f1 = 0
prec = 0
rec = 0
clf = RandomForestClassifier(n_estimators=500, max_depth=50, random_state=1)
for i, (train, test) in enumerate(splits):
    probas_ = classifier.fit(x[train], y[train]).predict_proba(x[test])
    clf.fit(x[train], y[train])
    newp = clf.predict(x[test])
    predict_values = classifier.fit(x[train], y[train]).predict(x[test])
    score += accuracy_score(y[test], predict_values)
    scorenew += accuracy_score(y[test], newp)
    mccscore += mcc(y[test], predict_values)
    # Compute ROC curve and area the curve

    # print ("f1score:", f1_score(y[test], predict_values))
    print(classification_report(y[test], predict_values))
    f1 += f1_score(y[test], predict_values)
    fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
    mean_tpr += interp(mean_fpr, fpr, tpr)
    mean_tpr[0] = 0.0
    prec = prec + precision_score(y[test], predict_values)
    rec = rec + recall_score(y[test], predict_values)
    roc_auc = auc(fpr, tpr)
    print("f1score:", f1_score(
        y[test], predict_values), 2 * ((prec * rec) / (prec + rec)))
    pl.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))

pl.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
print("score:", score / 5, "mcc:", mccscore / 5,
      "precision:", prec / 5, "recall:", rec / 5)
mean_tpr /= len(splits)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
print('Mean AUC:', mean_auc)
print("Mean f1:", f1 / 5)
print("RF", scorenew / 5)
print(clf.feature_importances_)

pl.plot(mean_fpr, mean_tpr, 'k--',
        label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

pl.xlim([-0.05, 1.05])
pl.ylim([-0.05, 1.05])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('LBTope')
# pl.show()

if args.test != 0:  # to test on blind data
    if args.test != '0': #condition to check if the testing data is given or not
        newdata = open(args.test, 'r')
        anew = []
        for l in newdata.readlines():
            if l[0] == '#':
                continue
            else:
                anew.append(l.strip().split())
        for i in range(len(anew)):
            aapscore = aap(anew[i][0], aapdic, args.average, args.reduce)
            aatscore = aat(anew[i][0], aatdic, args.average, args.reduce)
            anew[i] = np.insert(anew[i], len(anew[i]) - 1, aapscore)
            anew[i] = np.insert(anew[i], len(anew[i]) - 1, aatscore)
            '''aac = AAC(anew[i][0])
            aac = [AAC(anew[i][0])[x]
                   for x in ['A', 'F', 'G', 'I', 'L', 'N', 'P', 'Q', 'T', 'W']]
            # print(aac)
            anew[i] = np.insert(anew[i], len(anew[i]) - 1, aac)
            # anew[i] = np.insert(anew[i], len(anew[i]) - 1, aac)

            dpc = DPC(anew[i][0])
            anew[i] = np.insert(anew[i], len(anew[i]) - 1, dpc)'''
        anew = np.array(anew)
        #print(anew)
        ynewtest = anew[:, -1].astype(int)
        xnewtest = anew[:, 1:-1].astype(float)
        xnewtest = np.array(preprocessing.scale(xnewtest))
        probas_ = classifier.fit(x, y).predict_proba(xnewtest)

        predict_values = classifier.fit(x, y).predict(xnewtest)
        score = accuracy_score(ynewtest, predict_values)
        scorenew=accuracy_score(ynewtest,clf.fit(x,y).predict(xnewtest))
        mccscore = mcc(ynewtest, predict_values)
        print(classification_report(ynewtest, predict_values))
        print("score:", score, "mcc:", mccscore, "newscore with rf:", scorenew)
        print("auc:", roc_auc_score(ynewtest, probas_[:, 1]))
        print("test:", auc)
        testp = open('predictions.txt','w')
        for i in range(len(anew)):
            print(anew[i, 0] + "\t" + str(anew[i, -1]) + "\t" + str(predict_values[i]) + "\t" + str(probas_[i]), end='\n', file=testp)
        testp.close()


