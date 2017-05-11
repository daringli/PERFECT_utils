from profile import Profile

class Generator(object):
    
    def __repr__(self):
        d = vars(self)
        l = []
        keys=d.keys()
        for k in keys:
            l.append(k + " : " + str(d[k]))
        return "\n".join(l)
    

    def __str__(self):
        return str(vars(self))

    def generate():
        # example of a minimal implementation
        p = Profile(lambda: None)
        p.generator_dict = vars(self)
    
