#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 16:05:11 2021

@author: comelassarat
"""

import math
import sys
import time

start_time = time.time()

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


########################################## END DATA ######################################################

# Scaling algorithm from "A Revealing Introduction to Hidden Markov Models by Mark Stamp"

def A_and_B(matrix_A, matrix_B, matrix_pi, emissions):
    
   



    
    C = [0]*nb_of_emissions
    the_alphas=[0]*nb_of_emissions
    c = 0
    alpha=[0]*N
    emission_0 = int(emissions[0])
    for i in range(N):
        alpha[i] = matrix_pi[0][i]*matrix_B[i][emission_0]
        c = c + alpha[i]

    c = 1/c
    C[0] = c
    for i in range(N):
        alpha[i] = alpha[i]*c
    #print(alpha)
    the_alphas[0] = alpha
    
    
    for t in range(1,nb_of_emissions):
        alpha=[0]*N
        c=0
        emission_t = int(emissions[t])
        #print(len(emissions))
        #print(alpha)
        for i in range(N):
            alpha_i = 0
            for j in range(N):
                alpha_i = alpha_i + the_alphas[t-1][j]*matrix_A[j][i]
            alpha_i = alpha_i*matrix_B[i][emission_t]
            alpha[i] = alpha_i
            c = c + alpha_i
        c = 1/c
        C[t] = c
        for i in range(N):
            alpha[i] = alpha[i]*c
        the_alphas[t] = alpha
    
    

    ############# BETA ################
    
    


    beta_T_moins_1 = [0]*N
    beta = [0]*nb_of_emissions
    for i in range(N):
        beta_T_moins_1[i] = C[-1]
    beta[-1] = beta_T_moins_1

    beta_t_plus_1 = beta_T_moins_1
    for t in reversed(range(nb_of_emissions-1)):
        emission_t = int(emissions[t+1])
        beta_t= [0]*N
        for i in range(N):
            beta_t_i = 0
            for j in range(N):
                beta_t_i = beta_t_i + matrix_A[i][j]*matrix_B[j][emission_t]*beta_t_plus_1[j]
            beta_t_i = beta_t_i*C[t]
            beta_t[i] = beta_t_i
        #print(beta_t)
        beta[t] = beta_t
        beta_t_plus_1 = beta_t
    #print(beta)
    #beta.reverse()


######################## GAMMA ############################
    
    
    di_Gamma_function = {}
    Gamma_function=[0]*nb_of_emissions
    

    for t in range(0,nb_of_emissions-1):  #We go to T-2
        #print(t)
        gamma=[0]*N
        emission_t= int(emissions[t+1])
        #denom = 0
        for i in range(N):
            gamma_t_i = 0
            #numer = 0
            alpha_t_i = the_alphas[t][i]
            for j in range(N):
                gamma_ij = alpha_t_i*matrix_A[i][j]*matrix_B[j][emission_t]*beta[t+1][j]
                gamma_t_i = gamma_t_i + gamma_ij
                di_Gamma_function[i,j,t] = gamma_ij 
                #numer+= gamma_ij
            gamma[i]=gamma_t_i
            #denom += gamma_t_i

            
        Gamma_function[t] =gamma
    Gamma_function[nb_of_emissions-1] = the_alphas[-1]

    
    
########################## ESTIMATION ##############################


    
    #Estimation of pi
    matrix_pi = [Gamma_function[0]]
    
    #Re-estimate of matrix_A
    for i in range(N):
        denom_A = 0
        for t in range(nb_of_emissions-1):
            denom_A = denom_A + Gamma_function[t][i]
        for j in range(N):
            numer_A = 0
            for t in range(nb_of_emissions-1):
                numer_A = numer_A + di_Gamma_function[i,j,t]
            matrix_A[i][j] = numer_A/denom_A
    
    
    #Re-estimate of matrix_B
    for i in range(N):
        denom_B = 0
        for t in range(nb_of_emissions-1):
            denom_B = denom_B + Gamma_function[t][i]
        for j in range(K):
            numer_B = 0
            for t in range(nb_of_emissions-1):
                if emissions[t] == j:
                    numer_B = numer_B + Gamma_function[t][i]
            matrix_B[i][j] = numer_B/denom_B


    logProb = 0
    for i in range(nb_of_emissions):
        logProb = logProb + math.log(C[i])
    logProb = -logProb
    

    return matrix_A, matrix_B, matrix_pi, logProb


maxIters = 1000
oldLogProb = -float('inf')
iters=1

matrix_A, matrix_B, matrix_pi, logProb = A_and_B(matrix_A_extracted, matrix_B_extracted, matrix_pi_extracted, emissions_extracted)

while (iters < maxIters and logProb > oldLogProb and time.time() - start_time < 0.9):
    oldLogProb = logProb
    matrix_A, matrix_B, matrix_pi, logProb = A_and_B(matrix_A, matrix_B, matrix_pi, emissions_extracted)
    iters+=1

matrix_A_good_format = [len(matrix_A),len(matrix_A)]
matrix_B_good_format = [len(matrix_B),len(matrix_B[0])]

for i in range(len(matrix_A)):
    matrix_A_good_format = matrix_A_good_format + matrix_A[i]   
    
for i in range(len(matrix_B)):
    matrix_B_good_format = matrix_B_good_format + matrix_B[i]   
    

print(' '.join(str(l) for l in matrix_A_good_format))
print(' '.join(str(l) for l in matrix_B_good_format))





