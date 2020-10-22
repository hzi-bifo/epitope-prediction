import sys
import score_fasta
import argparse
import os.path
import warnings
warnings.filterwarnings("ignore")

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
        parser.add_argument("-m", "--model", default = "general", help=(" general for general predictions and viral for viral predictions"))
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


    aap_file=''
    aat_file=''
    print("model is:",model)
    if model.strip() == "general":
        print("Making predictions using general model")
        aap_file="aap-general.normal"
        aat_file="aat-general.normal"
        model = 'svm-general.pickle'
    if model.strip() == "viral":
        print("Making predictions using viral model")
        aap_file="aap-viral.normal"
        aat_file="aat-viral.normal"
        model = 'svm_viral.pickle'

    peptide_list, pred_probability = score_fasta.scoremodel(ifile, "./model/"+model, aap_file, aat_file )
    out = open("./output/"+outputfile, 'w')
    epitopelist=[]
    nonepitopelist=[]
    
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            epitopelist.append(i)
        '''if pred_probability[i][1] < 0.5:
            nonepitopelist.append(i)
            print("Peptides predicted as non-epitopes:")
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]))
            print("Non-epitopes \n", file=out)
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]), file=out)'''
    
    if len(epitopelist) > 0:
        print("Peptides predicted as epitopes:")
        print("Epitopes \n", file=out)
        for i in epitopelist:
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]))
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]), file=out)
    else:
        print("No peptides were predicted as epitopes")
    out.close()

