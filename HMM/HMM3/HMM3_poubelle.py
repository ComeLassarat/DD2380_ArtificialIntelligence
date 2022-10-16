#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:03:32 2021

@author: comelassarat
"""

'''
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
            
 ''' 
 
       

'''
def list_of_alphas_and_betas(matrix_A, matrix_B, matrix_pi, emissions):
#def list_of_alphas_and_betas(data):
    #matrix_A, matrix_B, matrix_pi, emissions = decoding_from_text(data)
    nb_of_emissions = len(emissions)
    C = []
    #print(nb_of_emissions)
    #print(matrix_B)


   ############################ ALPHAS ##############################
    
    #Initialization
    
    b = [[matrix_B[i][int(emissions[0])] for i in range(len(matrix_B))]] #extracting the column of B matching witg first emission
    alpha = product_of_lines(b,matrix_pi) #first alpha
    
    
    c0 = 1/sum(alpha[0])        #because alpha = [[x,y,z,...]]
    C.append(c0)
    alpha_scaled = []
    for i in range(len(alpha[0])):
        alpha_scaled_i = alpha[0][i]*c0
        alpha_scaled.append(alpha_scaled_i)
    
    the_alphas = [alpha_scaled]

    
    #Bucle
    for i in range(1,nb_of_emissions):
        emission = int(emissions[i])
        c = 0
        b = [[matrix_B[j][emission] for j in range(len(matrix_B))]]
        t_A_alpha = matrix_prod([alpha_scaled],matrix_A)
        alpha = product_of_lines(t_A_alpha,b)

        #Scaling
        c = 1/sum(alpha[0])
        #print('alpha')
        #print(alpha)
        #print(sum(alpha[0])) 
        C.append(c)
        alpha_scaled = []
        for i in range(len(alpha[0])):
            alpha_scaled_i = alpha[0][i]*c
            alpha_scaled.append(alpha_scaled_i)

        the_alphas.append(alpha_scaled) 
    
     ########################## BETAS ############################
    

    N = len(matrix_A)   #number of states
    emissions.reverse()
    the_betas = []
    
    #Initialisation
    
    beta_t_plus_1 = [[C[-1] for i in range(N)]]
    the_betas.append(beta_t_plus_1[0])
    
    #Bucle
    
    for em in emissions[0:nb_of_emissions-1]:          #reversed emissions
        #print('em')
        #print(em)
        count = nb_of_emissions - 2 #To compute c_t
        new_beta = []       
        b = [[matrix_B[i][int(em)] for i in range(len(matrix_B))]]
        #print('b')
        #print(b)
        
        
        for i in range(N):
            sum_ = 0
            for j in range(N):
                element = beta_t_plus_1[0][j]*b[0][j]*matrix_A[i][j]
                sum_ = sum_ + element
                
            new_beta.append(sum_)
        #print('new_beta')
        #print(new_beta)
        scaled_new_beta = [new_beta[n]*C[count] for n in range(len(new_beta))] #c_t computed
        #print('scaled_new_beta')
        #print(scaled_new_beta)
            
        count = count -1
        the_betas.append(scaled_new_beta)
        beta_t_plus_1 = [scaled_new_beta]
        
    #print(len(the_betas))
    
    ################### RETURN RESULTS ######################
    
    return the_alphas, the_betas, C
'''

 '''for t in range(0,nb_of_emissions-1):  #We go to T-2
        #print(t)
        #b = [[matrix_B[i][int(emissions[nb_of_emissions - n-1-1])] for i in range(len(matrix_B))]]  # we go from emissions T-1 to 1
        b = [[matrix_B[i][int(emissions[t+1])] for i in range(len(matrix_B))]]
        #sum_alpha_T = 0
        gamma = []
        #for k in range(N):
            #sum_alpha_T = sum_alpha_T + the_alphas[-1][k]
            #print(the_alphas[-1][k])
            
        
        for i in range(N):
            sum_ = 0
            for j in range(N):
                #print(the_alphas[int(emissions[t])][i])
                #print(the_betas[int(emissions[t+1])][j])
                gamma_ij = the_alphas[int(emissions[t])][i]*matrix_A[i][j]*b[0][j]*the_betas[int(emissions[t+1])][j]#/sum_alpha_T
                di_Gamma_function[i,j,t] = gamma_ij
                sum_ = sum_ + gamma_ij
            
            gamma.append(sum_)   #gamma vector at emissions t
            
        Gamma_function.append(gamma)    #All the gamma vectors
        
        
        
    Gamma_function.append(the_alphas[-1])
    #print(len(Gamma_function))
    
    
    
    return Gamma_function , di_Gamma_function
                
'''  


'''

def a(Gamma_function, di_Gamma_function, matrix_A, nb_of_emissions):
    
    N = len(matrix_A)
    
    for i in range(N):
        sum_=0
        for s in range(nb_of_emissions-1):
            sum_ = sum_ + Gamma_function[s][i]
        
        for j in range(N):
            sum_di = 0
            for n in range(nb_of_emissions-1):
                sum_di = sum_di +  di_Gamma_function[i,j,n]
            matrix_A[i][j] = sum_di/sum_
    
    return matrix_A
                

def b(Gamma_function, di_Gamma_function, matrix_B, emissions):
    N = len(matrix_B)
    K = len(matrix_B[0])
    nb_of_emissions = len(emissions)
    
    
    for j in range(N):
        sum_=0
        for s in range(nb_of_emissions-1):
            sum_ = sum_ + Gamma_function[s][j]
        
        for k in range(K):
            sum_2 = 0
            for n in range(nb_of_emissions-1):
                if k == emissions[n]:
                    sum_2 = sum_2 +  Gamma_function[n][j]
                matrix_B[j][k] = sum_2/sum_
                
        return matrix_B
                
    
'''
    
