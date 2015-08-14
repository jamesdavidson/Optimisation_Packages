#This is the main program used to coordinate all the other
import VectCalc, Test_Functs, Optimisers
import os
import numpy as np
import itertools as itt

upper_dimension_limit = 6
univ_low_limit = -2
univ_up_limit = 2.2
lim_inc = 0.5

data_directory = r'/home/senex/Documents/Graduate Diploma/Operations Research/Group Assignment/Data'

class function_Block(object):
    '''
    classdocs
    '''
    def __init__(self,Funct): 
        '''Function properties object,Grad properties object, Hessian Properties Object'''
        self.Funct = Funct #This object you put in x and out comes the function value
        Funct.Gen_Vars()
    def expression_give(self):
        return self.Funct.Gen_Funct()
    def variables_give(self):
        return self.Funct.Gen_Vars()
    def grad_give(self):
        Grd = VectCalc.Grad(self.Funct.Gen_Funct(),self.Funct.Gen_Vars()) #This object you put in x and out comes the grad value
        return Grd
    def hessian_give(self):
        Hessian = VectCalc.Hessian(self.grad_give(), self.Funct.Gen_Vars()) #This object you put in x and out comes the hessian value
        return Hessian

for dmsn in range(2,upper_dimension_limit+1): #All sets of n-th dimensional Rosenbrocks
    opt_Ros = Test_Functs.Rosenbrock(dmsn) #Creates a Rosenbrock object
    Ros_Comb = function_Block(opt_Ros)
    for x in itt.product(np.arange(univ_low_limit,univ_up_limit,lim_inc),repeat=dmsn):
        Data = open(os.path.join(data_directory,'%d-th dimensional Newton Rosenbrock data.txt' %(dmsn)),'a+')
        x0=np.array(x)
        [x_min, iter_no, time] = Optimisers.Newton(x0,Ros_Comb)
        Data.write(str(list(x0)) + ';' + str(list(x_min)) + ';' + str(iter_no) + ';' + str(time) + '\r\n')
        Data.close()
