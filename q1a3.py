from math import exp,log
from MNL import *
#Monte Carlo function
def _using_MC(s,f,N,a,m):
    val = 0
    x = np.array(mu_li_co_ge(s, a, m, N))/m
    for i in range(N):
        val = val + f(x[i])
    ans = val/N
    return ans

def p(x,a):
    return a * exp(-x)

#Defining the integration function
def int(x):
    return exp(-x**2)

#Defining the function for importnace sampling
#Given distribution after integration from 0 to 1 = 1.
#this gives $\alpha$ = e/e-1 and sample(x) = -ln(1-x/alpha).
def isam(x):
    return int(-log(1-x))/p(x,1)


#Using Monte carlo of integration without importance sampling
#Giving parameters for MLCG
# s = 1, a = 572 , m = 16381, N = 10000
val1 = _using_MC(1,int,10000,572,16381)
val2 = _using_MC(1,isam,10000,572,16381)

#----Printing Results----#
print("Value of Integral without importance sampling : ",val1)
print("Value of Integral with importance sampling : ",val2)

#Value of Integral without importance sampling :  0.7465264360450006
#Value of Integral with importance sampling :  0.7552730221984496