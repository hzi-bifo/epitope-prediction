import sys
import score_fasta
import argparse
import os
import warnings
warnings.filterwarnings("ignore")


def error(er):
    sys.exit(er)

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message+"\n")
        self.print_help()
        sys.exit(2)

    def parse_args(self, args=None, namespace=None):
        namespace = super().parse_args(args, namespace)
        if namespace.length < 4:
            MyParser.error(self, "Peptide length can't be less than 4. Please enter a value >= 4. For best results use 20(default).")
        if namespace.cutoff > 1:
            MyParser.error(self, "Cut-off score can't be more than 1. Please enter a value between 0 and 1.")
        if namespace.cutoff < 0:
            MyParser.error(self, "Cut-off score can't be less than 0. Please enter a value between 0 and 1.")
        return namespace


def options():
    
    parser = MyParser(
    description=("Takes an input fasta file and outputs the score of the peptides "),
    formatter_class=argparse.HelpFormatter,
    usage=("%(prog)s [-h]"
           " -i inputfile -o outputfolder [-m modeltype] [-l length] [-c cutoff"), add_help=False)
    group_required = parser.add_argument_group('required arguments')
    group_required.add_argument("-i", "--input",  help=(" input fasta file.Please put the file in input folder or enter the full path \n"), required=True)
    group_required.add_argument("-o", "--output", help=(" output folder where predictions will be saved \n"), required=True)
    group_optional = parser.add_argument_group('optional arguments')
    group_optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='show this help message and exit')
    group_optional.add_argument("-m", "--model", default = "general",choices=['general', 'viral'], help=(" general for general predictions(default) and viral for viral predictions\n"))
    group_optional.add_argument("-l", "--length", default = 20, type=int, help=("length of the peptides. Should be an integer >=4. For best results, use peptide length = 20(default).\n"),metavar= '>=4')
    group_optional.add_argument("-c", "--cutoff", default = 0.5, type=float, help=("cut-off score (default = 0.5) for a peptide to be considered an epitope. Should be a number between 0 and 1.\n"), metavar = '0.0 to 1.0')
    '''group_optional.add_argument("-c", "--cutoff", default = 0.5, type=float, help=("cut-off score (default = 0.5) for a peptide to be considered an epitope. Should be a number between 0 and 1.\n"), metavar = '0.0 to 1.0')'''
    #args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    args = parser.parse_args()
    #print(args.input)
    inputfile = args.input
    outputfile = args.output
    model = args.model
    length = args.length
    cutoff = args.cutoff
    return inputfile, outputfile, model, length, cutoff

if __name__ == "__main__":
    inputfile, outputdir, model, length, cutoff = options()
    
    #print("model is:",model)
    if model.strip() == "general":
        print("\nPredictions will be made using general model \n")
        model = 'svm-ibce.pickle'
    if model.strip() == "viral":
        print("\nPredictions will be made using viral model \n")
        model = 'svm-viral.pickle'

    if os.path.isfile('./protvec/sp_sequences_4mers_vec.bin') == False:
        error("Error. Protvec binary file is missing in the protvec folder. See README file")

    ifile=''
    if os.path.isfile('./input/'+inputfile) == True: # check if the inputfile is in the input folder
        print("Reading inputfile: /input/"+inputfile+"\n")
        ifile = "./input/"+inputfile
    elif os.path.isfile(inputfile)==True: #check if the inputfile exists in the path
        print ("Reading inputfile: "+inputfile+"\n")
        ifile = inputfile
    else:
        error("Error. Inputfile not found in the specified location. Please enter the correct name/path of the inputfile.")

    ofolder=''
    outputdir = os.path.realpath(outputdir)
    if os.path.isdir(outputdir) == True:
        print("Output folder already exists. Predictions will be put in the",os.path.basename(outputdir),"folder. \n")
        ofolder = outputdir
    if os.path.isdir(outputdir) == False and os.path.isdir(os.path.dirname(outputdir)) == True :
        print("Output folder:",os.path.basename(outputdir)," will be created at the specified path:",os.path.dirname(outputdir),"\n")
        os.mkdir(outputdir)
        ofolder = outputdir
    if os.path.isdir(outputdir) == False and os.path.isdir(os.path.dirname(outputdir)) == False :
        error("Error. Unable to find or create the output folder. Please check the output path")

    sequences, pep, peploc, seqid, pred = score_fasta.scoremodel(ifile, "./model/"+model, length )

    
    score_fasta.print_file(sequences, pep, peploc, seqid, pred, cutoff, length, ofolder)

    if len(sequences.keys())==1:
        score_fasta.print_stdout(pep, peploc, seqid, pred, cutoff, length)
    else:
        print("Fasta file has more than one sequence. The predictions won't be displayed here. Please check the predictions in the \'"+ofolder+"\' folder. \n")
    print("Done!")
