

def nplusone(a,c,m,n):    
    return (a*n + c)%m              

a = 1664525
b= 2**32
c,n = 1013904223, 0        

for i in range(10):    
        n = nplusone(a,c,b,n)
        print n
