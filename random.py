
def nplusone(a=1664525,c=1013904223,m=2**32,n=0):
    return (a*n+c)%m
n0=0
for i in range(10):    
    n1 = nplusone(a=1664525,c=1013904223,m=2**32,n=n0)
    n0=n1         
    print n1
