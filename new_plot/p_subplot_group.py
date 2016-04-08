from p_subplot import perfect_subplot
import numpy

class perfect_subplot_group:
    def __init__(self,p_subplot_list,groups,get_all=False,logic="or"):
        if get_all is False:
            if logic=="or":
                self.p_subplot_list = [x for x in p_subplot_list if not set(x.groups).isdisjoint(groups)]
            if logic=="and":
                self.p_subplot_list = [x for x in p_subplot_list if set(groups).issubset(x.groups)] 
        if get_all is True:
            self.p_subplot_list = p_subplot_list

    def setattrs(self,attr_name,value):
        #sets an attribute of all members in the list
        for p_subplot in self.p_subplot_list:
            setattr(p_subplot,attr_name,value)
        
    def getattrs(self,attr_name):
        #returns a list of the attribute for all members in the list
        return [getattr(p_subplot,attr_name) for p_subplot in self.p_subplot_list]


    def get_max(self,attr_name):
        maxes = [numpy.max(x) for x in self.getattrs(attr_name)]
        return numpy.max(maxes)

    def get_min(self,attr_name):
        mines = [numpy.min(x) for x in self.getattrs(attr_name)]
        return numpy.min(mines)

    def set_middle_attr(self,attr_name,value):
        #sets an attribute of the middle element of a list.
        #This is designed for plots where subplots in the same group are close to each other
        #so that it makes sense to only label one of them.
        l=len(self.p_subplot_list)
        i=int(l/2.0)
        setattr(self.p_subplot_list[i],attr_name,value)
