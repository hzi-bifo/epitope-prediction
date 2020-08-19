import sys
import score_fasta
import argparse
import os.path

def options():
    try:
        parser = argparse.ArgumentParser(
        description=("Takes an input fasta file and outputs the score of peptides "
                     "additonal parameters for testing model"),
        formatter_class=argparse.RawTextHelpFormatter,
        usage=("%(prog)s inputfile [-h]"
               "[-i inputfile][-m modelfile][-o outputfile]"))
        parser.add_argument("-i", "--input",  help=(" input fasta file"))
        parser.add_argument("-o", "--output", help=(" output file with predicitons"))
        parser.add_argument("-m", "--model", default = "svm-model.pickle", help=(" trained machine learning model"))
        args = parser.parse_args()

        inputfile = args.input
        outputfile = args.output
        model = args.model
    except:
        print("Invalid parameters. See help")
        sys.exit()
    return inputfile, outputfile, model

if __name__ == "__main__":
    inputfile, outputfile, model = options()
    ifile=''
    if os.path.isfile('./input/'+inputfile)==True: # check if the inputfile is in the input folder
        print("Reading inputfile: /input/"+inputfile)
        ifile = "./input/"+inputfile
    elif os.path.isfile(inputfile)==True: #check if the inputfile exists in the path
        print ("Reading inputfile: "+"inputfile")
        ifile = inputfile
    else:
        print ("Error in reading inputfile. Please enter the full path of the inputfile.")
        sys.exit()

    peptide_list, pred_probability = score_fasta.scoremodel(ifile, "./model/"+model )
    out = open("./output/"+outputfile, 'w')
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]), file=out)
    out.close()

