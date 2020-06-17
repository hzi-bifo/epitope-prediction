# EpitopeVec: Linear Epitope Prediction Using DeepProtein Sequence Embeddings
EpitopeVec predicts linear b-cell epitopes. It takes a protein sequence (FASTA format) as input and then lists the peptides (of 20 amino acids length) that can be epitopes along with their respective predicition probability. It is based on a SVM model trained on a large set of experimentally verified epitopes and makes use of different amino acid features, k-mer features, and distributed vector reprsentation of proteins (ProtVec).  

## Requirements
* Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/hcv-mds 
  ```
* ```Python 3``` with the following packages:
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
   
## Usage
 Run the main.py file with the following command.
    
    python3 main.py -i inputfile.fasta -o outputfilename -m machine-learning-model
    
   

