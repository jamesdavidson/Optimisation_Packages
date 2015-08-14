'''
Created on 23/05/2014

@author: senex
'''

def var_matcher(varbs,vals):
    match_List = []
    for i in range(0,len(varbs)):
        match_List.append((varbs['%d' %(i+1)],vals[i]))
    return match_List

def multiVariableHalfOpen(x0, expr, varbs, dk, T=0.01):
    '''x0 is the starting coordinates, function is an object with an extractable expression object (expression)
    and an extractable variable object (vbl), dk is the descent direction vector and T is the increment parameter'''
    k=1
    p = x0
    q = x0 + T*dk
    fp = expr.subs(var_matcher(varbs,p))
    fq = expr.subs(var_matcher(varbs,q))
    while fp > fq:
        k += 1
        p=q
        fp=fq
        q=p+(2**(k-1))*T*dk
        fq = expr.subs(var_matcher(varbs,q))
    if k == 1:
        a = 0
        b = T
    elif k == 2:
        a = 0
        b = 3*T
    else:
        u = range(0,k)
        v = range(0,k-2)
        for i in range(0,len(v)):
            v[i]=2**v[i]
        for i in range(0,len(u)):
            u[i]=2**u[i]
        a = T*sum(v)
        b = T*sum(u)
    return [a,b]
