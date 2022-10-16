#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 10:22:34 2021

@author: comelassarat
"""

import sys

def matrix_op(Y,Z):  #A&B list matrix --> A*B
        n = len(Y) #lines
        m = len(Y[0]) #rows
        p = len(Z) #lines
        l = len(Z[0]) #rows

        if m != p :
            return False
        
        prod = []

        for i in range (n):
            line = []
            for j in range (l):
                sum = 0
                for k in range (m):
                    sum += Y[i][k] * Z[k][j]
                line.append(sum)
            prod.append(line)
        return prod
    
''''
def decoding_from_text(data):
    A=[]
    B=[]
    pi=[]
    lambdaa = [A,B,pi]
    h = open(data, 'r')
    content = h.readlines()
    nb_line=-1
    for line in content:
        nb_line+=1
        for i in line:
            lambdaa[nb_line] = line.split()
            lambdaa[nb_line] = [float(lambdaa[nb_line][k]) for k in range(len(lambdaa[nb_line]))]



    A = lambdaa[0][2:]
    N = int(lambdaa[0][0])    #number of possible states
    B = lambdaa[1][2:]
    #print('B')
    #print(B)
    K = int(lambdaa[1][1])    #number of possible observations
    pi = lambdaa[2][2:] #initial probability
    

    matrix_A = []
    matrix_B = []
    matrix_pi = [pi]
    #print('matrix_pi')
    #print(matrix_pi)
    

    #matrix

    for i in range (N):
        matrix_A.append(A[i*N:(i+1)*N])
        matrix_B.append(B[i*K:(i+1)*K])
    #print('matrix_B')
    #print(matrix_B)
    
    return matrix_A,matrix_B,matrix_pi
'''

A_str = sys.stdin.readline()
B_str = sys.stdin.readline()
pi_str = sys.stdin.readline()

A=[]
B=[]
C=[]

lambdaa = [A,B,C]

matrix_str = [A_str,B_str,pi_str]
nb_line=-1
for line in matrix_str:
    nb_line+=1
    for i in line:
        lambdaa[nb_line] = line.split()
        lambdaa[nb_line] = [float(lambdaa[nb_line][k]) for k in range(len(lambdaa[nb_line]))]


A = lambdaa[0][2:]
N = int(lambdaa[0][0])    #number of possible states
B = lambdaa[1][2:]
#print('B')
#print(B)
K = int(lambdaa[1][1])    #number of possible observations
pi = lambdaa[2][2:] #initial probability


matrix_A = []
matrix_B = []
matrix_pi = [pi]
#print('matrix_pi')
#print(matrix_pi)


#matrix

for i in range (N):
    matrix_A.append(A[i*N:(i+1)*N])
    matrix_B.append(B[i*K:(i+1)*K])

def HMM0():
    #matrix_A,matrix_B,matrix_pi = decoding_from_text(data)
    pi_A = matrix_op(matrix_pi,matrix_A) #pi*A
    result = matrix_op(pi_A,matrix_B) #pi*A*B
    result[0].insert(0,len(result[0])) #adding the dimensions
    result[0].insert(0,len(result)) #adding the dimensions
    return ' '.join(str(l) for l in result[0])

print(HMM0())
    


    

