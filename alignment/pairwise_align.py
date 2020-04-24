#! /usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 17:12:34 2015
@authors: Kaveh Karbasi & Ehsaneddin Asgari
"""
import sys
sys.path.append('../')
import os.path
import timeit
import re


class PWAlignment:
    def __init__(self, seq1, seq2, seq_type):
        '''
            the class constructor gets:
            (1) the first and (2)the second sequence and (3)the scoring matrix
        '''
        # Characters score matrix , gap penalty , input characters
        # this method sets the self.scores, self.characters,self.gap_pnlty of this class
        if seq_type == 'dna':
            self.create_scores_dict("../config/dna_dist.txt")
        elif seq_type == 'rna':
            self.create_scores_dict("../config/rna_dist.txt")
        elif seq_type == 'protein':
            self.create_scores_dict("../config/protein_dist.txt")
        else:
            print("Not a valid type")
            exit()
        # Sequences
        self.seq1 = seq1.upper()
        self.seq2 = seq2.upper()

        # check validity
        if not (self.checkSeqValidity(self.seq1) or self.checkSeqValidity(self.seq2)):
            print("Characters in the sequence does not match the scoring matrix")
            exit()

        # Alignment score matrix
        self.matrix = [[[0, 0] for i in self.seq1 + '1'] for j in self.seq2 + '1'];
        # Calculate the alignmet score matrix
        self.calc_matrix();

    def create_scores_dict(self, mat_file_path):
        '''
        Creates a dictionary of scores of all pair characters
        returns the pair character score dictionary along with gap penalty and 
        allowable characters
        '''
        # open the score matrix file
        infile = open(mat_file_path)
        matData = infile.read();
        # init the score matrix valuse
        self.scores = {}
        lines = matData.split("\n")
        # detemining the characters in the scoring matrix
        self.characters = lines[0].split();

        # ******* Error handing
        # check if the header does have any numerical character
        if (any(i.isdigit() for i in lines[0])):
            print("Incorrect format detected in the scoring matrix:\n ** no numerical character is allowed")
            exit()

        N = len(self.characters);
        # ******* Error handing
        # check if the number of lines is consistent with the number of characters
        if not len(lines) == N + 2:
            print("Incorrect format detected in the scoring matrix :\n ** # of lines doesn't  match the character set")
            exit()
        # setting the score matrix values
        for lindex, line in enumerate(lines[1:-1]):
            try:
                vals = [int(x) for x in line.split()];
            except:
                print("Incorrect format detected in the scoring matrix in line: " + str(
                    lindex + 2) + ":\n ** only numerical value is allowed")
                exit()
            # ******* Error handing
            # detecting the inconsistency between # of char and the int in each row of matrix
            if not (len(vals) == N):
                print("Incorrect format detected in the scoring matrix in line: " + str(lindex + 2))
                exit()
            for cindex, char in enumerate(self.characters):
                self.scores["".join([self.characters[lindex], char])] = vals[cindex]

        # parsing the gap penalty
        # ******* Error handing
        try:
            vals = [int(x) for x in lines[-1].split()];
        except:
            print("Incorrect format detected in the scoring matrix: \n ** incorrect gap penalty: a single number ")
            exit()
        # ******* Error handing
        # if more than one value is specified
        if not (len(vals) == 1):
            print("Incorrect format detected in the scoring matrix: \n ** incorrect gap penalty: a single number")
            exit()
        self.gap_pnlty = int(vals[0])

    def checkSeqValidity(self, sequence):
        '''
            This method checks if the sequence read from the fasta file
            matches the characters specified in the scoring matrix
        '''
        # compares if the characters in the sequence are are subset of chars in characters
        if set(sequence) <= set(self.characters):
            return True;
        else:
            print(set(sequence))
            return False;

    def calc_matrix(self):
        '''
        Calculates the alignmet score matrix
        dynamic programming
        '''
        J = range(len(self.matrix[0]))
        I = range(len(self.matrix))

        # initialization of the first column and the first row
        for i in I:
            self.matrix[i][0][0] = self.gap_pnlty * i
            self.matrix[i][0][1] = 2
        for j in J:
            self.matrix[0][j][0] = self.gap_pnlty * j
            self.matrix[0][j][1] = 1

        # following the dynamic programming rules for pairwise alignment
        for i in I[1:]:
            for j in J[1:]:
                # we need to calculate three options and calculate the optimum to assign to the current cell
                Match = self.matrix[i - 1][j - 1][0] + self.scores["".join([self.seq1[j - 1], self.seq2[i - 1]])]
                Insert = self.matrix[i][j - 1][0] + self.gap_pnlty
                Delete = self.matrix[i - 1][j][0] + self.gap_pnlty
                # 0 is diagonal, 1 is horizantal, 2 is vertical
                pathList = [Match, Insert, Delete]
                # assign the best value
                self.matrix[i][j][0] = max(pathList)
                # keep the pointer to the previous cell
                self.matrix[i][j][1] = pathList.index(self.matrix[i][j][0]);

    def getScoreMatrix(self):
        '''
        Alignmet score matrix getter
        '''
        return self.matrix;

    def getAlignScore(self):
        '''
        Alignmet score matrix getter
        '''
        return self.matrix[-1][-1][0];

    def printAlignment(self):
        '''
           This function iteratively go from the buttom of the table
           to the head to find the aligned sequences and finally
           print the aligned sequences in chunks of 80 characters
        '''
        J = len(self.matrix[0]) - 1
        I = len(self.matrix) - 1
        # aligned version of the first sequence
        s1 = ""
        # aligned version of the second sequence
        s2 = ""
        # in this loop we start from the final score and track the
        # path to the starting point and incrementally build the aligned sequences
        # in s1 and s2

        while not (I == 0 and J == 0):
            # find the previous cell
            alignPath = self.matrix[I][J][1];
            # Previous cell is the diagonal cell
            if alignPath == 0:
                # no gap penalty
                # update the aligned sequences
                s1 = self.seq1[J - 1] + s1
                s2 = self.seq2[I - 1] + s2
                # update the pointer to the current cell
                I = I - 1
                J = J - 1
            # Previous cell is on the left
            elif alignPath == 1:
                # update the aligned sequences
                s1 = self.seq1[J - 1] + s1
                s2 = '-' + s2
                # update the pointer to the current cell
                J = J - 1
            # previous cell is on the top of the current
            else:
                # update the aligned sequences
                s1 = '-' + s1
                s2 = self.seq2[I - 1] + s2
                # update the pointer to the current cell
                I = I - 1
            # now we need to print them in chunks of 80
        alignment_output = ""
        s1_chunks = [s1[i:i + 80] for i in range(0, len(s1), 80)]
        s2_chunks = [s2[i:i + 80] for i in range(0, len(s1), 80)]
        chunk_num = 0;
        # preparing the alignment output by adding pair of alignment in the length of 80
        for s1_chunk in s1_chunks:
            alignment_output += s1_chunk + "\n"
            alignment_output += s2_chunks[chunk_num] + "\n\n"
            chunk_num = chunk_num + 1
        # We print everything once to be efficient in time
        print(alignment_output.rstrip())

    # Print the scoring matrix
    # it prints every lines together to save time
    def printScoreMatrix(self):
        row_num = 0;
        # first row
        matrix_out = "|*|*|" + ('|'.join(self.seq1)) + "|\n"
        # the second sequence in on the vertical header
        ver_header = "*" + self.seq2;
        for row in self.matrix:
            # getting the scorse for the current row
            s = [i[0] for i in row]
            # adding the fist colums
            s.insert(0, '|' + ver_header[row_num])
            row_num += 1
            # join othe columns
            matrix_out += ('|'.join([str(i) for i in s])) + "|\n"
        # print the scoring matrix
        print(matrix_out)
