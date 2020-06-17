import sys
import score_fasta
import argparse

def options():
    try:
        parser = argparse.ArgumentParser(
        description=("Takes an input fasta file and outputs the score of peptides "
                     "additonal parameters for testing model"),
        formatter_class=argparse.RawTextHelpFormatter,
        usage=("%(prog)s inputfile [-h]"
               "[-i inputfile][-m modelfile][-o outputfile]"))
        parser.add_argument("-i", "--input", help=(" input fasta file"))
        parser.add_argument("-o", "--output", help=(" output file with predicitons"))
        parser.add_argument("-m", "--model", help=(" trained machine learning model"))
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
    peptide_list, pred_probability = score_fasta.scoremodel(inputfile, model )
    out = open("./output/"+outputfile, 'w')
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]), file=out)
    out.close()

