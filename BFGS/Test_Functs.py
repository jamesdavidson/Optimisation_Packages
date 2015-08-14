'''
Created on 19/05/2014

@author: senex
'''
import sympy
from math import *

class Rosenbrock(object):
    '''
    classdocs
    '''

    def __init__(self, dims): #Creates instantiated dimensions and empty variable dictionary
        '''
        Constructor
        '''
        self.dims=dims
        self.vbl = {}
        
    def Gen_Vars(self): #Creates list from x1 to x_n
        for i in range(1,self.dims+1):
            self.vbl['%d'%(i)]= sympy.symbols('x%d' %(i))
        return self.vbl
            
    def Gen_Funct(self): #Generates the n-th dimensional Rosenbrock function as a symbolic expression
        self.expression = 0
        for i in range(1,self.dims):
            self.expression = self.expression \
                                + 100*((self.vbl['%d' %(i+1)]-(self.vbl['%d' %(i)])**2)**2)\
                                +(1-self.vbl['%d' %(i)])**2
        return self.expression

class QuadSimp(object):
    '''
    classdocs
    '''

    def __init__(self, dims): #Creates instantiated dimensions and empty variable dictionary
        '''
        Constructor
        '''
        self.dims=dims
        self.vbl = {}
        
    def Gen_Vars(self): #Creates list from x1 to x_n
        for i in range(1,self.dims+1):
            self.vbl['%d'%(i)]= sympy.symbols('x%d' %(i))
        return self.vbl
            
    def Gen_Funct(self): #Generates the n-th dimensional Rosenbrock function as a symbolic expression
        self.expression = 0
        for i in range(1,self.dims+1):
            self.expression = self.expression + (3.5%i)*(self.vbl['%d' %(i)]-i%2)**2
        return self.expression

class Styb(object):
    '''
    Based upon the Styblinski-Tang Function
    '''

    def __init__(self, dims): #Creates instantiated dimensions and empty variable dictionary
        '''
        Constructor
        '''
        self.dims=dims
        self.vbl = {}
        
    def Gen_Vars(self): #Creates list from x1 to x_n
        for i in range(1,self.dims+1):
            self.vbl['%d'%(i)]= sympy.symbols('x%d' %(i))
        return self.vbl
            
    def Gen_Funct(self): #Generates the n-th dimensional Rosenbrock function as a symbolic expression
        self.expression = 0 #cos(self.vbl['x1' %(i)])
        for i in range(1,self.dims+1):
            self.expression = self.expression + (self.vbl['%d' %(i)])**4 - 16*(self.vbl['%d' %(i)])**2 + 5*(self.vbl['%d' %(i)]) 
        return self.expression
