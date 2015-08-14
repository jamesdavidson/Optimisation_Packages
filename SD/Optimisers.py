'''
Created on 23/05/2014

@author: senex
'''
import Bound_Find, OneD_Search
import time
import numpy as np
from math import sqrt
import sys

def var_matcher(varbs,vals):
    match_List = []
    for i in range(0,len(varbs)):
        match_List.append((varbs['%d' %(i+1)],vals[i]))
    return match_List

def normAr(lister):
    tally = 0
    for element in lister:
        tally = tally + element**2
    norm = sqrt(tally)
    return norm

def clear_all_vars():
    sys.modules[__name__].__dict__.clear()

def Steepest_Desc(x0, function_Block, iter_limit=1000, threshold=0.01):
    '''Function block object must have a Grad in object and a function
    value calculator given a coordinate. Hessian can be empty.
    This works for symbolic functions and black box models'''
    Expr = function_Block.expression_give()
    Varbs = function_Block.variables_give()
    Grad = function_Block.grad_give()

    pointsSD_file = open(r'/home/senex/Documents/2D Steepest Descent Point list.txt','a+')
    start_Time = time.time()
    xk = x0
    print xk
    k=0
    dk = function_Block.grad_give()
    var_Val_List = var_matcher(Varbs,xk)
    for i in range(0,len(Grad)):
        dk[i]=-Grad[i].subs(var_Val_List)
    Grad_Mag = normAr(dk)
    pointsSD_file.write('\r\n'+str(xk)+'\r\n')
    while Grad_Mag>threshold and k < iter_limit:
        k+=1
        [a,b] = Bound_Find.multiVariableHalfOpen(xk, Expr, Varbs, dk)
        step = OneD_Search.multiVariableGoldenSectionSearch(Expr, Varbs,a,b,xk,dk)
        xk = xk + step*dk
        print xk
        pointsSD_file.write(str(xk)+'\r\n')
        Grad = function_Block.grad_give()
        var_Val_List = var_matcher(Varbs,xk)
        dk = function_Block.grad_give()
        for i in range(0,len(dk)):
            dk[i]=-Grad[i].subs(var_Val_List)
        Grad_Mag = normAr(dk)

    x_min = xk

    if Grad_Mag>threshold and k==iter_limit:
        iter_no = None
        time_taken = None
    else:
        iter_no = k
        time_taken = time.time() - start_Time

    pointsSD_file.close()

    return [x_min,iter_no,time_taken]

def Newton(x0, function_Block, iter_limit=1000, threshold=0.01):
    Expr = function_Block.expression_give()
    Varbs = function_Block.variables_give()
    Grad = function_Block.grad_give()

    pointsN_file = open(r'/home/senex/Documents/2D Newton Point list.txt','a+')
    start_Time = time.time()
    k=0
    xk = x0
    print xk

    var_Val_List = var_matcher(Varbs,xk)
    for i in range(0,len(Grad)):
        Grad[i]=-Grad[i].subs(var_Val_List)
    Grad_Mag = normAr(Grad)

    Hessian = function_Block.hessian_give()
    for i in range(0,Hessian.shape[0]):
        for j in range(0,Hessian.shape[1]):
            Hessian[i,j] = float(Hessian[i,j].subs(var_Val_List))
    Hessian = np.array(Hessian,dtype=np.float32)/float(np.amax(Hessian))

    pointsN_file.write('\r\n'+str(xk)+'\r\n')
    while Grad_Mag>threshold and k < iter_limit:
        k += 1
        if (np.linalg.eig(Hessian)[0]>np.zeros(Hessian.shape[1])).all():
            dk = np.dot(-np.linalg.inv(Hessian),Grad)
        else:
            dk = -Grad
        [a,b] = Bound_Find.multiVariableHalfOpen(xk, Expr, Varbs, dk)
        step = OneD_Search.multiVariableGoldenSectionSearch(Expr, Varbs,a,b,xk,dk)
        xk = xk + step*dk
        print xk
        pointsN_file.write(str(xk)+'\r\n')
        Grad = function_Block.grad_give()
        var_Val_List = var_matcher(Varbs,xk)
        for i in range(0,len(Grad)):
            Grad[i]=Grad[i].subs(var_Val_List)
        Grad_Mag = normAr(Grad)
        Hessian = function_Block.hessian_give()
        for i in range(0,Hessian.shape[0]):
            for j in range(0,Hessian.shape[1]):
                Hessian[i,j] = Hessian[i,j].subs(var_Val_List)
        Hessian = np.array(Hessian,dtype=np.float32)/float(np.amax(Hessian))

    x_min = xk

    if Grad_Mag>threshold and k==iter_limit:
        iter_no = None
        time_taken = None
    else:
        iter_no = k
        time_taken = time.time() - start_Time

    pointsN_file.close()

    return [x_min,iter_no,time_taken]

def BFGS(x0, function_Block, H0=None, iter_limit=1000, threshold=0.01):
    Expr = function_Block.expression_give()
    Varbs = function_Block.variables_give()
    Grad = function_Block.grad_give()
    if H0 == None:
        H0=np.eye(len(Grad),dtype=np.float32)

    pointsQ_file = open(r'/home/senex/Documents/2D BFGS Point list.txt','a+')
    start_Time = time.time()
    k=0
    xk = x0
    print xk
    xk_old = x0
    H_old = H0

    var_Val_List = var_matcher(Varbs,xk)
    for i in range(0,len(Grad)):
        Grad[i]=Grad[i].subs(var_Val_List)
    Grad_Mag = normAr(Grad)

    pointsQ_file.write('\r\n'+str(xk)+'\r\n')
    while Grad_Mag>threshold and k < iter_limit:
        k += 1
        H0 = np.array(H0,dtype=np.float32)/float(np.amax(H0))
        dk = -np.dot(H_old,Grad)
        [a,b] = Bound_Find.multiVariableHalfOpen(xk, Expr, Varbs, dk)
        step = OneD_Search.multiVariableGoldenSectionSearch(Expr, Varbs,a,b,xk,dk)
        xk = xk + step*dk
        print xk
        pointsQ_file.write(str(xk)+'\r\n')
        xk_new = xk_old + step*dk
        sk = xk_new - xk_old

        Grad_New = function_Block.grad_give()
        new_var_Val_List = var_matcher(Varbs,xk_new)
        for i in range(0,len(Grad_New)):
            Grad_New[i]=Grad_New[i].subs(new_var_Val_List)

        Grad_Old = function_Block.grad_give()
        old_var_Val_List = var_matcher(Varbs,xk_old)
        for i in range(0,len(Grad_Old)):
            Grad_Old[i]=Grad_Old[i].subs(old_var_Val_List)

        gk = Grad_New - Grad_Old
        rk = np.dot(H_old,gk)/np.dot(sk,gk)
        H_new = H_old + ((1+np.dot(rk,gk))/np.dot(sk,gk))*np.outer(sk,sk) - np.outer(sk,rk) - np.outer(rk,sk)

        xk_old = xk_new
        H_old = H_new

        Grad_Mag = normAr(Grad_New)

    x_min = xk

    if Grad_Mag>threshold and k==iter_limit:
        iter_no = None
        time_taken = None
    else:
        iter_no = k
        time_taken = time.time() - start_Time

    pointsQ_file.close()

    return [x_min,iter_no,time_taken]
