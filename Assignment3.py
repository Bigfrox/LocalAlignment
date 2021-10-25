'''
Bio Computing Assignment 3, Local Pairwise Sequence Alignment
2016253072
명수환(Myeong Suhwan)
없는 문자 : BJOUXZ
'''

import re
import random
from datetime import datetime
import string
import sys
import numpy as np

def getDataFromFile(filename):
    data1 = ""
    data2 = ""
    try:
        with open(filename, 'r') as file:
            line = None
            line = file.readline() # comment
            if not line: # * input_file is empty
                print("[-] No protein sequence . . .")
                exit(1)
            
            if(line[0] == '>'): # * must accept only first FASTA format.
                print("1st FASTA : [Comment]")
                print(line)
            else:
                print("[-] No correct format . . .")
                exit(1)
            while line != '':
                line = file.readline()
                if line:
                    
                    if(line[0] == '>'): # * accept 2nd FASTA format.
                        print("2nd FASTA : [Comment]")
                        print(line)
                        while line != '':
                            
                            line = file.readline()
                            
                            if line:
                                if(line[0] == '>'): # * ignore 3rd and after FASTA format.
                                    print("[+] 3rd FASTA is ignored . . .")
                                    return data1, data2
                                else:
                                    data2 += line.strip('\n')

                    else:
                        
                        data1 += line.strip('\n')
                else:
                    print("[-] Need one more sequence.")
                    exit(1)
                
    except FileNotFoundError:
        print("[-] No input file . . .")
        
    return data1, data2


def processing(data):
    data = data.replace(" ", "")
    data = data.replace("\t", "")
    data = data.upper()

    return data

def GetPenalty(alphabet1, alphabet2):
    matrix = [9, -1, -1, -3, 0, -3, -3, -3, -4, -3, -3, -3, -3, -1, -1, -1, -1, -2, -2, -2,-1, 4, 1, -1, 1, 0, 1, 0, 0, 0, -1, -1, 0, -1, -2, -2, -2, -2, -2, -3,-1, 1, 4, -1, 0, -2, 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, 0, -2, -2, -2,-3, -1, -1, 7, -1, -2, -2, -1, -1, -1, -2, -2, -1, -2, -3, -3, -2, -4, -3, -4, 0, 1, 0, -1, 4, 0, -2, -2, -1, -1, -2, -1, -1, -1, -1, -1, 0, -2, -2, -3, -3, 0, -2, -2, 0, 6, 0, -1, -2, -2, -2, -2, -2, -3, -4, -4, -3, -3, -3, -2, -3, 1, 0, -2, -2, 0, 6, 1, 0, 0, 1, 0, 0, -2, -3, -3, -3, -3, -2, -4, -3, 0, -1, -1, -2, -1, 1, 6, 2, 0, -1, -2, -1, -3, -3, -4, -3, -3, -3, -4, -4, 0, -1, -1, -1, -2, 0, 2, 5, 2, 0, 0, 1, -2, -3, -3, -2, -3, -2, -3,-3, 0, -1, -1, -1, -2, 0, 0, 2, 5, 0, 1, 1, 0, -3, -2, -2, -3, -1, -2, -3, -1, -2, -2, -2, -2, 1, -1, 0, 0, 8, 0, -1, -2, -3, -3, -3, -1, 2, -2, -3, -1, -1, -2, -1, -2, 0, -2, 0, 1, 0, 5, 2, -1, -3, -2, -3, -3, -2, -3, -3, 0, -1, -1, -1, -2, 0, -1, 1, 1, -1, 2, 5, -1, -3, -2, -2, -3, -2, -3, -1, -1, -1, -2, -1, -3, -2, -3, -2, 0, -2, -1, -1, 5, 1, 2, 1, 0, -1, -1, -1, -2, -1, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3, 1, 4, 2, 3, 0, -1, -3, -1, -2, -1, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2, 2, 2, 4, 1, 0, -1, -2, -1, -2, 0, -2, 0, -3, -3, -3, -2, -2, -3, -3, -2, 1, 3, 1, 4, -1, -1, -3, -2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 6, 3, 1, -2, -2, -2, -3, -2, -3, -2, -3, -2, -1, 2, -2, -2, -1, -1, -1, -1, 3, 7, 2, -2, -3, -2, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3, 1, 2, 11]

    matrix = np.reshape(matrix,(20,20))

    blosum = ['C','S','T','P','A','G','N','D','E','Q','H','R','K','M','I','L','V','F','Y','W']
    
    for i in range(len(blosum)):
        if alphabet1 == blosum[i]:
            alphabet1 = i
            break
            

    for j in range(len(blosum)):
        if alphabet2 == blosum[j]:
            alphabet2 = j
            break
            
    
    
    penalty = matrix[i][j]
    return penalty

def LocalAlignmentUsingDP(data1, data2):
    
    #! 'i' is index for data 1
    #! 'j' is index for data 2
    m = len(data1)+1
    n = len(data2)+1
    max_score = 0
    max_i,max_j = 0,0
    
    score = [[0 for col in range(n)] for row in range(m)]
    
    for i in range(m):
        score[i][0] = 0 # * fill these blocks with 0
    
    for j in range(n):
        score[0][j] = 0 # * fill these blocks with 0
    for i in range(1,m):
        for j in range(1,n):
            
            score[i][j] = max(score[i-1][j]+gap_penalty,
            score[i][j-1]+gap_penalty,
            score[i-1][j-1]+GetPenalty(data1[i-1],data2[j-1]))
            if max_score < score[i][j]:
                max_score = score[i][j]
                max_i = i
                max_j = j
                
                
    #!DEBUG
    # for i in range(m):
    #     print(score[i])

    subsequence1 = []
    subsequence2 = []
    count = score[m-1][n-1] # * length of LCS
    
    print("[*] Local Alignment Score : ", max_score)
    #print("position : ",max_i,max_j)
    m -=1
    n -=1
    
    #* print("[*] Start the Backtracking . . .")
    subsequence1,subsequence2 = BackTracking(data1,data2,score,max_i,max_j,subsequence1,subsequence2)
    #* print("[*] End of the Backtracking . . .")
    return subsequence1,subsequence2

def BackTracking(data1,data2,score, row, col,subsequence1,subsequence2):
    
    m = row
    n = col

    if m <= 0 or n <= 0:
        return subsequence1,subsequence2
    if score[m][n] == 0:
        return subsequence1,subsequence2

    #? max score means where to go
    #? have to add '-' when moved horizontally or vertically.
    direction = max(score[m][n-1]+gap_penalty, score[m-1][n]+gap_penalty, score[m-1][n-1]+GetPenalty(data1[m-1],data2[n-1]))
    if direction == score[m-1][n-1]+GetPenalty(data1[m-1],data2[n-1]): #* diagonal
        subsequence1.insert(0,data1[m-1])
        subsequence2.insert(0,data2[n-1])
        m -=1
        n -=1
        subsequence1,subsequence2 = BackTracking(data1,data2,score,m,n,subsequence1,subsequence2)
    elif direction == score[m][n-1]+gap_penalty: #* horizontal - add gap '-' to data1 side
        subsequence1.insert(0,'-') #* data1 side
        subsequence2.insert(0,data2[n-1])
        n -= 1
        subsequence1,subsequence2 = BackTracking(data1,data2,score,m,n,subsequence1,subsequence2)

    elif direction == score[m-1][n]+gap_penalty: #* vertical - add gap '-' to data2 side
        subsequence1.insert(0,data1[m-1])
        subsequence2.insert(0,'-') #* data2 side
        m -= 1
        subsequence1,subsequence2 = BackTracking(data1,data2,score,m,n,subsequence1,subsequence2)

    return subsequence1,subsequence2

def main():
    if len(sys.argv) != 2:
        print("No input file.")
        print("<Usage> assignment3.py input_filename.txt")
        return -1;
    
    input_filename = sys.argv[1]
    
    global gap_penalty
    gap_penalty = -5
    #* 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    data_not = re.compile(r'[^acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY]')
    data1 = ""
    data2 = ""
    
    data1,data2 = getDataFromFile(input_filename)
    data1 = processing(data1)
    data2 = processing(data2)

    if not data1 and not data2: # * if input data is empty
        print("No Protein Sequences . . .")
        exit(1)
    if not data1 or not data2:
        print("[!] Need one more sequence / No correct format ")
        exit(1)

    no_protein = data_not.findall(data1+data2)
    if no_protein:
        print("No Protein Sequences . . .")
        print(no_protein)
        exit(1)

    print("Data1 : ",data1)
    print("Data2 : ",data2)
    print("\n")
    start_time = datetime.now()
    result1,result2 = LocalAlignmentUsingDP(data1,data2)
    result1 = "".join(result1)
    result2 = "".join(result2)
    
    index = 0
    while index < len(result1):
        print(result1[index:index+60])
        print(result2[index:index+60])
        index += 60
        print("\n")
    
    
    print("\n")
    

    #print(datetime.now())
    print("[+] Time Elapsed : ", datetime.now() - start_time, "microseconds")

if __name__ == '__main__':
    main()