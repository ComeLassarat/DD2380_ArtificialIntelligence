#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 22:25:48 2021

@author: comelassarat
"""

import sys


def matrix_prod(A,B):  #A&B list matrix --> A*B
        n = len(A) #lines
        m = len(A[0]) #rows
        p = len(B) #lines
        l = len(B[0]) #rows

        if m != p :
            return False
        
        prod = []

        for i in range (n):
            line = []
            for j in range (l):
                sum = 0
                for k in range (m):
                    sum += A[i][k] * B[k][j]
                line.append(sum)
            prod.append(line)
        return prod


def product_of_lines(U,V):   #U and V are vectors as following : [[a,b,c,d,...]]
    l = len(U[0]) #or len(V[0]) because they have the same length
    prod = []
    for i in range(l):
        prod.append(U[0][i]*V[0][i])

    return [prod]

def transpose(A):
    t_A = []
    nb_of_lines = len(A)
    nb_of_columns = len(A[0])
    for i in range(nb_of_columns):
            new_line = [A[j][i] for j in range(nb_of_lines)]
            t_A.append(new_line)
    return t_A

########################################## DATA ######################################################


A_str = sys.stdin.readline()
B_str = sys.stdin.readline()
pi_str = sys.stdin.readline()
emissions_str = sys.stdin.readline()

A=[]
B=[]
pi=[]
emissions= []

lambdaa = [A,B,pi,emissions]

matrix_str = [A_str,B_str,pi_str,emissions_str]
nb_line=-1
for line in matrix_str:
    nb_line+=1
    for i in line:
        lambdaa[nb_line] = line.split()
        lambdaa[nb_line] = [float(lambdaa[nb_line][k]) for k in range(len(lambdaa[nb_line]))]


A = lambdaa[0][2:]
N = int(lambdaa[0][0])    #number of possible states

B = lambdaa[1][2:]
K = int(lambdaa[1][1])    #number of possible observations

pi = lambdaa[2][2:] #initial probability

emissions = lambdaa[3][1:]


matrix_A = []
matrix_B = []
matrix_pi = [pi]
#print('matrix_pi')
#print(matrix_pi)


#matrix

for i in range (N):
    matrix_A.append(A[i*N:(i+1)*N])
    matrix_B.append(B[i*K:(i+1)*K])



##################################### END DATA ###########################################

def HMM1():       
    nb_of_emissions = len(emissions)

    
    #initialization
    
    b = [[matrix_B[i][int(emissions[0])] for i in range(len(matrix_B))]] #extracting the column of B matching witg first emission
    alpha = product_of_lines(b,matrix_pi) #alpha1
    
    #Bucle
    for i in range(1,nb_of_emissions):
        emission = int(emissions[i])
        b = [[matrix_B[j][emission] for j in range(len(matrix_B))]]
        t_A_alpha = matrix_prod(alpha,matrix_A)
        alpha = product_of_lines(t_A_alpha,b)
    
    return sum(alpha[0])
    

print(HMM1())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
