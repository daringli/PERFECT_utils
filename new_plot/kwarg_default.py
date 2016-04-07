def kwarg_default(key,default,**kwargs):
    if key in kwargs.keys():
        return kwargs[key]
    else:
        return default
            
if __name__=="__main__":
    d={'t':"non-default"}
    t=kwarg_default('t',"default",**d)
    print t
    s=kwarg_default('s',"default",**d)
    print s
