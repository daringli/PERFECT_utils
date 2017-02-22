from __future__ import division

class Profile(object):
    def __init__(self,f):
        self.f = f
    def __call__(self,x):
        return self.f(x)

    def __neg__(self): #unary minus
        return Profile(lambda x: -self.f(x))
    def __pos__(self): #unary plus
        return self

    def __add__(self,other):
         if callable(other):
             return Profile(lambda x: self(x) + other(x))
         else:
             return Profile(lambda x: self(x) + other)
    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other):
        return self.__add__(-other)
    def __rsub__(self,other):
        return self.__add__(-other)
    
    def __mul__(self,other):
        if callable(other):
            return Profile(lambda x: self(x) * other(x))
        else:
            return Profile(lambda x: self(x) * other)
    def __rmul__(self,other):
        return self.__mul__(other)
    
    def __truediv__(self,other):
        if callable(other):
            return Profile(lambda x: self(x) / other(x))
        else:
            return Profile(lambda x: self(x) / other)

    def __rtruediv__(self,other):
        if callable(other):
            return Profile(lambda x: other(x) / self(x))
        else:
            return Profile(lambda x: other / self(x))


if __name__ == "__main__":
    import numpy
    x=numpy.linspace(0,2,3)
    
    # operator tests between Profiles
    p1 = Profile(lambda x : x)
    p2 = Profile(lambda x : x**2)
    p3 = p1 + p2 # x + x^2
    p4 = p1 * p2 # x^3
    p5 = p1 / p2 # 1/x
    p6= -p1

    print p1(x) # 0 1 2
    print p2(x) # 0 1 4
    print p3(x) # 0 2 6
    print p4(x) # 0 1 8
    print p5(x) # nan 1 0.5
    print p6(x) # 0 -1 -2

    print "==========================="
    # operator tests between Profiles and functions
    p1 = Profile(lambda x : x)
    p2 = lambda x : x**2
    p3 = p1 + p2 # x + x^2
    p4 = p1 * p2 # x^3
    p5 = p1 / p2 # 1/x
    p6= -p1

    print p1(x) # 0 1 2
    print p2(x) # 0 1 4
    print p3(x) # 0 2 6
    print p4(x) # 0 1 8
    print p5(x) # nan 1 0.5
    print p6(x) # 0 -1 -2

    print "==========================="
    # operator tests between functions and Profiles
    p1 = lambda x : x
    p2 = Profile(lambda x : x**2)
    p3 = p1 + p2 # x + x^2
    p4 = p1 * p2 # x^3
    p5 = p1 / p2 # 1/x
    p6= lambda x : -x

    print p1(x) # 0 1 2
    print p2(x) # 0 1 4
    print p3(x) # 0 2 6
    print p4(x) # 0 1 8
    print p5(x) # nan 1 0.5
    print p6(x) # -0 -1 -2

    print "==========================="
    # operator tests between numbers and Profiles
    p1 = 1.2
    p2 = Profile(lambda x : x**2)
    p3 = p1 + p2 # x + x^2
    p4 = p1 * p2 # x^3
    p5 = p1 / p2 # 1/x

    print p1 # 1.2
    print p2(x) # 0 1 4
    print p3(x) # 1.2 2.2 5.2
    print p4(x) # 0 1.2 4.8
    print p5(x) # inf 1.2 0.3

    print "==========================="
    # operator tests between Profiles and numbers
    p1 = Profile(lambda x : x)
    p2 = 1.2 
    p3 = p1 + p2 # x + x^2
    p4 = p1 * p2 # x^3
    p5 = p1 / p2 # 1/x

    print p1(x) # 0 1 2
    print p2 # 1.2
    print p3(x) # 1.2 2.2 3.2
    print p4(x) # 0 1.2 2.4
    print p5(x) # 0 0.83333 1.666666

