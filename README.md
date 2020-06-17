# EpitopeVec: Linear Epitope Prediction Using DeepProtein Sequence Embeddings
EpitopeVec predicts linear b-cell epitopes. It takes a protein sequence (FASTA format) as input and then lists the peptides (of 20 amino acids length) that can be epitopes along with their respective predicition probability. It is based on a SVM model trained on a large set of experimentally verified epitopes and makes use of different amino acid features, k-mer features, and distributed vector reprsentation of proteins (ProtVec).  

## Requirements
* Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/hcv-mds 
  ```
* **```Python 3```** with the following packages:
    * **numpy 1.16.4**
    * **scipy 1.4.1**
    * **matplotlib 3.1.3**
    * **sklearn 0.22.1**
    * **pydpi 1.0**
    * **biopython 1.71.0**
    
   
  If these are not installed, you can install them with ``` pip ```. 
    ```
   pip3 install -r requirements.txt
   ```
   
 * Binary file for ProtVec representation of proteins. It can be downloaded from  and then it should be put in the ```protvec``` folder.
   
## Usage
 Run the main.py file with the following command.
    
    python3 main.py -i inputfile.fasta -o outputfilename -m machine-learning-model
    
 -i : This takes the input file with the protein sequence (in FASTA format) on which epitopes are to be predicted.
 
 -o : The name of the output file. This will contain the list of peptides which are predicted as epitopes.
  
 -m : Machine-learning model trained on the training set as a pickle file.
 
## Input
  The program takes a protein sequence as input. The input file should be in **FASTA** file format and should only use the 20 aa codes.       Please put the input file in the __input__ folder. Please take a look in the __input__ folder to see an example of input fasta       file: ```example.fasta``` 
  
## Output
The output file is created in the __output__ folder with the outputname provided. It is a txt file with two columns. First column is contains the amino acid sequences of the predicted peptides and second column is the predicted probability score for the likelihood of that peptide being an epitope.
    
   

