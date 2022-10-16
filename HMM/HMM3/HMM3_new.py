#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 16:05:11 2021

@author: comelassarat
"""

import math
import sys
import time

start_timestamp = time.time()

########################################## DATA ######################################################


A_str = sys.stdin.readline()
B_str = sys.stdin.readline()
pi_str = sys.stdin.readline()
emissions_str = sys.stdin.readline()

A = A_str.split(' ')
B = B_str.split(' ')
pi = pi_str.split(' ')
the_emissions = emissions_str.split(' ')

len_A = len(A)
len_B = len(B)
len_pi = len(pi)
len_emissions = len(the_emissions)

A=[float(A[i]) for i in range(len_A-1)]
B=[float(B[i]) for i in range(len_B-1)]
pi=[float(pi[i]) for i in range(len_pi-1)]
the_emissions=[float(the_emissions[i]) for i in range(len_emissions-1)]

N = int(A[0])
K = int(B[1])
nb_of_emissions = int(the_emissions[0])
#print(nb_of_emissions)

emissions_extracted = the_emissions[1:]
#print(len(emissions_extracted))
matrix_pi_extracted = [pi[2:]]

matrix_A_extracted = [0]*N
matrix_B_extracted = [0]*N

for i in range(N):
    matrix_A_extracted[i] = A[2+i*N:2+(i+1)*N]
    matrix_B_extracted[i] = B[2+i*K:2+(i+1)*K]

##################################### END DATA ###########################################











# Scaling algorithm from "A Revealing Introduction to Hidden Markov Models by Mark Stamp"

def list_of_alphas_and_betas2(matrix_A, matrix_B, matrix_pi, emissions):
#def list_of_alphas_and_betas2(data):
    
    #start1_timestamp = time.time()
    ############### ALPHA ###################
    #matrix_A, matrix_B, matrix_pi, emissions = decoding_from_text(data)
    
    #nb_of_emissions = len(emissions)
    C = [0]*nb_of_emissions
    the_alphas=[0]*nb_of_emissions
    c = 0
    alpha=[0]*N
    for i in range(N):
        alpha_i = matrix_pi[0][i]*matrix_B[i][int(emissions[0])]
        alpha[i] = alpha_i
        c = c + alpha_i

    c = 1/c
    C[0] = c
    for i in range(N):
        alpha[i] = alpha[i]*c
    #print(alpha)
    the_alphas[0] = alpha
    
    
    for t in range(1,nb_of_emissions):
        alpha=[0]*N
        c=0
        #print(alpha)
        for i in range(N):
            alpha_i = 0
            for j in range(N):
                alpha_i = alpha_i + the_alphas[t-1][j]*matrix_A[j][i]
                #print(alpha_i)
            alpha_i = alpha_i*matrix_B[i][int(emissions[t])]
            alpha[i] = alpha_i
            c = c + alpha_i
        c = 1/c
        C[t] = c
        for i in range(N):
            alpha[i] = alpha[i]*c
        the_alphas[t] = alpha
    
    ############# BETA ################
    beta_T_moins_1 = [0]*N
    beta = []
    for i in range(N):
        beta_T_moins_1[i] = C[-1]
    beta.append(beta_T_moins_1)
    
    for t in reversed(range(nb_of_emissions-1)):
        beta_t= [0]*N
        for i in range(N):
            beta_t_i = 0
            for j in range(N):
                beta_t_i = beta_t_i + matrix_A[i][j]*matrix_B[j][int(emissions[t+1])]*beta[0][j]
            beta_t_i = beta_t_i*C[t]
            beta_t[i] = beta_t_i
        beta.insert(0,beta_t)
    #print(len(C))
    #print(-start1_timestamp + time.time())
    return the_alphas, beta, C
        
  
def gamma(the_alphas, the_betas, matrix_A, matrix_B, matrix_pi, emissions):
#def gamma(data):
    #start1_timestamp = time.time()
    #print('###################################')
    
    #matrix_A, matrix_B, matrix_pi, emissions = decoding_from_text(data)
    
    
    #nb_of_emissions = len(emissions)
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
    print(nb_of_emissions)
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

    
    
    #nb_of_emissions = len(emissions)
    
    
    #Re-estimate of matrix_pi
    #new_matrix_pi = [[Gamma_function[0][i] for i in range(N)]]
    matrix_pi = [Gamma_function[0]]
    
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
        for j in range(K):
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

    

    #nb_of_emissions = len(emissions)
    maxIters = 4000
    
    oldLogProb = -1000000000

    
    the_alphas, the_betas, C = list_of_alphas_and_betas2(matrix_A, matrix_B, matrix_pi, the_emissions)
    Gamma_function , di_Gamma_function = gamma(the_alphas, the_betas, matrix_A, matrix_B, matrix_pi, the_emissions)
    matrix_A, matrix_B, matrix_pi = re_estimate(matrix_A, matrix_B, matrix_pi, Gamma_function, di_Gamma_function, the_emissions, C)
        
    logProb = 0
    for i in range(nb_of_emissions):
        logProb = logProb + math.log(C[i])
    logProb = -logProb
    
    iters = 1
    
    

    
    #print(start_timestamp - start1_timestamp)
    
    while (iters < maxIters and logProb > oldLogProb and time.time()- start_timestamp < 0.65):
        oldLogProb = logProb
        the_alphas, the_betas, C = list_of_alphas_and_betas2(matrix_A, matrix_B, matrix_pi, the_emissions)
        Gamma_function , di_Gamma_function = gamma(the_alphas, the_betas, matrix_A, matrix_B, matrix_pi, the_emissions)
        matrix_A, matrix_B, matrix_pi = re_estimate(matrix_A, matrix_B, matrix_pi, Gamma_function, di_Gamma_function, the_emissions, C)
            
        #print(iters)
        
        logProb = 0
        for i in range(nb_of_emissions):
            logProb = logProb + math.log(C[i])
        logProb = -logProb
        
        iters = iters + 1
        #print(time.time()- start_timestamp)
        #print(iters)
        #print('logProb')
        #print(logProb)
        #print('oldLogProb')
        #print(oldLogProb)


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
        
    

Baum_Welch(matrix_A_extracted, matrix_B_extracted, matrix_pi_extracted)
    
    


