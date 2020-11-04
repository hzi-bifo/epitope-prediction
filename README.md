# EpitopeVec: Linear Epitope Prediction Using DeepProtein Sequence Embeddings
---
## Contents
- [EpitopeVec Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

---


## EpitopeVec<a name="introduction"></a>
EpitopeVec predicts linear B-cell epitopes. It takes protein sequences (FASTA format) as input and then lists the peptides (of a given amino acids length) that can be epitopes along with their respective predicition probability. It is based on a SVM model trained on a large set of experimentally verified epitopes and makes use of different amino acid features, amino acid pair propensity scale, amino acid trimer propensity scale, and distributed vector representation of proteins (ProtVec).
The code for training new models and results from the EpitopVec article are available at https://github.com/hzi-bifo/epitope-prediction-paper  

## Requirements<a name="requirements"></a>

* **```Python 3```** with the following packages:
    * **numpy 1.18.5**
    * **scipy 1.5.0**
    * **matplotlib 3.2.2**
    * **scikit-learn 0.23.2**
    * **pydpi 1.0**
    * **biopython 1.77**
    * **tqdm 4.47.0**
    * **gensim 3.8.3**
    
   
  If these are not installed, you can install them with ``` pip ```. 
    ```
   pip3 install -r ./requirement/requirements.txt
   ```
   
  Additionally, **pydpi 1.0** from ```pip``` might be incompatible with **Python 3**. Please install the **pydpi** package from the provided ```requirement/pydpi.tar.gz``` file.
    ```
    pip3 install ./requirement/pydpi-1.0.tar.gz
    ```
   
 * Binary file for ProtVec representation of proteins can be downloaded using the following command in the ```protvec``` directory:
 
 ```
 cd protvec
 wget http://deepbio.info/embedding_repo/sp_sequences_4mers_vec.txt
 wget http://deepbio.info/embedding_repo/sp_sequences_4mers_vec.txt.bin -O sp_sequences_4mers_vec.bin
 ```
 
   
## Usage<a name="usage"></a>
* Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/epitope-prediction
  ```
* Run the main.py file with the following command.
  ```    
    python3 main.py -i inputfile.fasta -o outputfilename -m machine-learning-model (general OR viral)
  
   required arguments
    
    -i : This takes the input file with the protein sequences (in FASTA format) on which epitopes are to be predicted.
 
    -o : The name of the output directory. This will contain the list of peptides which are predicted as epitopes.
    
   optional arguments
     
     -m : Machine-learning model to use for making predictions. Use 'general' for general predictions(default) and 'viral' for viral prediction
     -l : Length of the peptides to be calculated from the sequence. Should be an integer >=4. For best results use length = 20(default)
     -c : cutoff score to consider for considering an peptide be epitope. Should be a number between 0 and 1. The default value is 0.5.
  ```
## Input<a name="input"></a>
  The program takes a file with protein sequences as input. The input file should be in **FASTA** file format and should only use the 20 aa codes. The fasta file can contain  a single or multiple protein sequences. If any sequence contains an invalid aa code or the length of the sequence is <4, then that sequence will be skipped while making predictions. Also, if the length of the sequence is < specified peptide length _(-l length)_, then the prediction will be made on the entire sequence.     
 Please put the input file in the _input_ folder or enter the full path of the inputfile. Please take a look in the _input_ folder to see an example of input fasta file: ```example.fasta``` 
  
## Output<a name="output"></a>
The output files are created in the _(-o output)_ folder. If the specified foler doesn't exist in the output path provided, it will be created. Files with the name of the sequence and suffix _.epitopes.txt_ for each sequence in the fasta file will be created in the output folder. The first column is the name of the sequence, second is the peptide, third is the starting postion of the peptide in the sequence, fourth is the ending postion of the peptide in the sequence and last column is the predicted probability score for the likelihood of that peptide being an epitope. By default, a predicted score of >= 0.5 is considered to be an epitope. To specify a different score cutoff, please use the _(-c cutoff)_ argument.
If the fasta file contains a single sequence then the predictions are also displayed on the commandline.

## License<a name="license"></a>
[GPL-3.0 License](https://github.com/hzi-bifo/epitope-prediction/blob/master/LICENSE.md)

## Citation<a name="citation"></a>
Add citation here.

## Contact<a name="contact"></a>
- Open an [issue](https://github.com/hzi-bifo/epitope-prediction/issues) in the repository.
- Send an email to Akash Bahai (akash.bahai@helmholtz-hzi.de).


