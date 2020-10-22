import sys
import score_fasta
import argparse
import os.path
import warnings
warnings.filterwarnings("ignore")

def options():
    try:
        parser = argparse.ArgumentParser(
        description=("Takes an input fasta file and outputs the score of peptides "),
        formatter_class=argparse.HelpFormatter,
        usage=("%(prog)s [-h]"
               "[-i inputfile][-o outputfile] [-m modeltype]"))
        parser.add_argument("-i", "--input",  help=(" input fasta file.Please put the file in input folder or enter full path \n"), required=True)
        parser.add_argument("-o", "--output", help=(" output file with predicitons \n"), required = True)
        parser.add_argument("-m", "--model", default = "general",choices=['general', 'viral'], help=(" general for general predictions and viral for viral predictions\n"))
        args = parser.parse_args()

        inputfile = args.input
        outputfile = args.output
        model = args.model
    except:
        sys.exit()
    return inputfile, outputfile, model

if __name__ == "__main__":
    inputfile, outputfile, model = options()
    
    aap_file=''
    aat_file=''
    #print("model is:",model)
    if model.strip() == "general":
        print("Predictions will be made using general model")
        aap_file="aap-general.normal"
        aat_file="aat-general.normal"
        model = 'svm-general.pickle'
    if model.strip() == "viral":
        print("Predictions will be made using viral model")
        aap_file="aap-viral.normal"
        aat_file="aat-viral.normal"
        model = 'svm-viral.pickle'

    ifile=''
    if os.path.isfile('./input/'+inputfile)==True: # check if the inputfile is in the input folder
        print("Reading inputfile: /input/"+inputfile)
        ifile = "./input/"+inputfile
    elif os.path.isfile(inputfile)==True: #check if the inputfile exists in the path
        print ("Reading inputfile: "+"inputfile")
        ifile = inputfile
    else:
        print ("Error in reading inputfile. Please enter the name/full path of the inputfile.")
        sys.exit()



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
        #print("Epitopes \n", file=out)
        for i in epitopelist:
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]))
            print(str(peptide_list[i])+"\t"+str(pred_probability[i][1]), file=out)
        print("The predicted epitopes are also available in the "+outputfile+" file in the output folder.")
    else:
        print("No peptides were predicted as epitopes")
    out.close()

