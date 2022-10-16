#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 16:05:11 2021

@author: comelassarat
"""

import math
import sys
import time

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











# Scaling algorithm from "A Revealing Introduction to Hidden Markov Models by Mark Stamp"

def list_of_alphas_and_betas2(matrix_A, matrix_B, matrix_pi, emissions):
#def list_of_alphas_and_betas2(data):
    
    #start1_timestamp = time.time()
    ############### ALPHA ###################
    #matrix_A, matrix_B, matrix_pi, emissions = decoding_from_text(data)
    N = len(matrix_A)
    nb_of_emissions = len(emissions)
    C = []
    the_alphas=[]
    c = 0
    alpha=[0]*N
    for i in range(N):
        alpha_i = matrix_pi[0][i]*matrix_B[i][int(emissions[0])]
        alpha[i] = alpha_i
        c = c + alpha_i

    c = 1/c
    C.append(c)
    for i in range(N):
        alpha[i] = alpha[i]*c
    #print(alpha)
    the_alphas.append(alpha)
    
    
    for t in range(1,nb_of_emissions):
        alpha=[]
        c=0
        #print(alpha)
        for i in range(N):
            alpha_i = 0
            for j in range(N):
                alpha_i = alpha_i + the_alphas[-1][j]*matrix_A[j][i]
                #print(alpha_i)
            alpha_i = alpha_i*matrix_B[i][int(emissions[t])]
            alpha.append(alpha_i)
            c = c + alpha_i
        c = 1/c
        C.append(c)
        for i in range(N):
            alpha[i] = alpha[i]*c
        the_alphas.append(alpha)
    
    ############# BETA ################
    beta_T_moins_1 = []
    beta = []
    for i in range(N):
        beta_T_moins_1.append(C[-1])
    beta.append(beta_T_moins_1)
    
    for t in reversed(range(nb_of_emissions-1)):
        beta_t= []
        for i in range(N):
            beta_t_i = 0
            for j in range(N):
                beta_t_i = beta_t_i + matrix_A[i][j]*matrix_B[j][int(emissions[t+1])]*beta[0][j]
            beta_t_i = beta_t_i*C[t]
            beta_t.append(beta_t_i)
        beta.insert(0,beta_t)
    #print(len(C))
    #print(-start1_timestamp + time.time())
    return the_alphas, beta, C
        
  
def gamma(the_alphas, the_betas, matrix_A, matrix_B, matrix_pi, emissions):
#def gamma(data):
    #start1_timestamp = time.time()
    #print('###################################')
    
    #matrix_A, matrix_B, matrix_pi, emissions = decoding_from_text(data)
    N = len(matrix_A)   #number of states
    
    nb_of_emissions = len(emissions)
    #emissions.reverse()
    
    #the_alphas, the_betas = list_of_alphas_and_betas(data)

    
    di_Gamma_function = {}
    #Gamma_function = [[0]*nb_of_emissions]
    Gamma_function=[0]*nb_of_emissions
    
 
   
    #print(the_betas)
    #print(the_alphas)
    #print(len(matrix_A[0]))
    #print(len(matrix_B[0]))
    #print('##########')
    #print(nb_of_emissions)
    #print(N)
    for t in range(0,nb_of_emissions-1):  #We go to T-2
        #print(t)
        gamma=[0]*N
        for i in range(N):
            gamma_t_i = 0
            for j in range(N):
                #print(i,j)
                gamma_ij = the_alphas[t][i]*matrix_A[i][j]*matrix_B[j][int(emissions[t+1])]*the_betas[t+1][j]
                gamma_t_i = gamma_t_i + gamma_ij
                di_Gamma_function[i,j,t] = gamma_ij 
            gamma[i]=gamma_t_i
            
        #gamma.append(the_alphas[-1])
        Gamma_function[t] =gamma
    Gamma_function[nb_of_emissions-1] = the_alphas[-1]
    

    
    return Gamma_function , di_Gamma_function
    



def re_estimate(matrix_A, matrix_B, matrix_pi, Gamma_function, di_gamma_function, emissions, C):
    
    start1_timestamp = time.time()

    N = len(matrix_A)
    M = len(matrix_B[0]) #number of observation symbols
    nb_of_emissions = len(emissions)
    
    
    #Re-estimate of matrix_pi
    #new_matrix_pi = [[Gamma_function[0][i] for i in range(N)]]
    #matrix_pi = [[Gamma_function[0][i] for i in range(N)]]
    
    #Re-estimate of matrix_A
    for i in range(N):
        denom = 0
        for t in range(nb_of_emissions-1):
            denom = denom + Gamma_function[t][i]
        for j in range(N):
            numer = 0
            for t in range(nb_of_emissions-1):
                numer = numer + di_gamma_function[i,j,t]
            matrix_A[i][j] = numer/denom
    
    
    #Re-estimate of matrix_B
    for i in range(N):
        denom = 0
        for t in range(nb_of_emissions-1):
            denom = denom + Gamma_function[t][i]
        for j in range(M):
            numer = 0
            for t in range(nb_of_emissions-1):
                if emissions[t] == j:
                    numer = numer + Gamma_function[t][i]
            matrix_B[i][j] = numer/denom
    
    #print(-start1_timestamp + time.time())
    return matrix_A, matrix_B, matrix_pi
    


def Baum_Welch(matrix_A, matrix_B, matrix_pi):       #data = file name
    #matrix_A,matrix_B,matrix_pi,emissions  = decoding_from_text(data)
    #start1_timestamp = time.time()

    start_timestamp = time.time()

    nb_of_emissions = len(emissions)
    maxIters = 100
    iters = 0
    oldLogProb = -float('inf')

    
    the_alphas, the_betas, C = list_of_alphas_and_betas2(matrix_A, matrix_B, matrix_pi, emissions)
    Gamma_function , di_Gamma_function = gamma(the_alphas, the_betas, matrix_A, matrix_B, matrix_pi, emissions)
    matrix_A, matrix_B, matrix_pi = re_estimate(matrix_A, matrix_B, matrix_pi, Gamma_function, di_Gamma_function, emissions, C)
        
    logProb = 0
    for i in range(nb_of_emissions):
        logProb = logProb + math.log(C[i])
    logProb = -logProb
    
    iters = iters + 1
    
    

    
    #print(start_timestamp - start1_timestamp)
        
    while (iters < maxIters and logProb > oldLogProb and time.time()- start_timestamp < 0.8):
        oldLogProb = logProb
        the_alphas, the_betas, C = list_of_alphas_and_betas2(matrix_A, matrix_B, matrix_pi, emissions)
        Gamma_function , di_Gamma_function = gamma(the_alphas, the_betas, matrix_A, matrix_B, matrix_pi, emissions)
        matrix_A, matrix_B, matrix_pi = re_estimate(matrix_A, matrix_B, matrix_pi, Gamma_function, di_Gamma_function, emissions, C)
            
        #print(iters)
        
        logProb = 0
        for i in range(nb_of_emissions):
            logProb = logProb + math.log(C[i])
        logProb = -logProb
        
        iters = iters + 1
        #print(time.time()- start_timestamp)
        #print(iters)
    

    matrix_A_good_format = [len(matrix_A),len(matrix_A)]
    matrix_B_good_format = [len(matrix_B),len(matrix_B[0])]
    
    for i in range(len(matrix_A)):
        matrix_A_good_format = matrix_A_good_format + matrix_A[i]   
        
    for i in range(len(matrix_B)):
        matrix_B_good_format = matrix_B_good_format + matrix_B[i]   
        

    matrix_A_str = ' '.join(str(l) for l in matrix_A_good_format)
    matrix_B_str = ' '.join(str(l) for l in matrix_B_good_format)

    print(matrix_A_str)
    print(matrix_B_str)
    #print(time.time()- start_timestamp)
        
    

Baum_Welch(matrix_A, matrix_B, matrix_pi)
    
    


