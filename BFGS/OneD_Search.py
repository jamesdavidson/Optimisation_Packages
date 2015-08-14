'''
Created on 23/05/2014

@author: senex
'''
from math import sqrt

def var_matcher(varbs,vals):
    match_List = []
    for i in range(0,len(varbs)):
        match_List.append((varbs['%d' %(i+1)],vals[i]))
    return match_List

def multiVariableGoldenSectionSearch(expr,varbs,a,b,xk,dk,tolerance=0.00001):
    gamma = (sqrt(5) - 1)/2
    k=1
    p = b - gamma*(b-a)
    q = a + gamma*(b-a)
    fp = expr.subs(var_matcher(varbs,xk + p*dk))
    fq = expr.subs(var_matcher(varbs,xk + q*dk))

    while (b-a) >= 2*tolerance:
        k += 1
        if fp <= fq:
            b = q
            q = p
            fq = fp
            p = b - gamma*(b-a)
            fp = expr.subs(var_matcher(varbs,xk+p*dk))
        else:
            a = p
            p = q
            fp = fq
            q = a + gamma*(b-a)
            fq=expr.subs(var_matcher(varbs,xk + q*dk))

    minEstimate = (a+b)/2;      #The calculated midpoint of the final interval
    #fminEstimate = expr.subs(var_matcher(varbs,xk + self.minEstimate*dk));
    return minEstimate
