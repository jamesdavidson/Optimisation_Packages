'''
Created on 21/05/2014

@author: senex
'''
from numpy import array
from sympy import diff

def Grad(expression, variables):
    pre_Sym_Grad = []
    for i in range(1,len(variables)+1):
        pre_Sym_Grad.append(diff(expression,variables['%d' %(i)]))
    Grd = array(pre_Sym_Grad)
    return Grd

def Hessian(grad_pre, variables):
    pre_Hessian = []
    for i in range(0,len(grad_pre)):
        pre_Hessian_Row = []
        grad_row = grad_pre[i]
        for j in range(1,len(variables)+1):
            pre_Hessian_Row.append(diff(grad_row,variables['%d'%(j)]))
        pre_Hessian.append(pre_Hessian_Row)
    Hess = array(pre_Hessian)
    return Hess
