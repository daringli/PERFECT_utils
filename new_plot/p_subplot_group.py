from p_subplot import perfect_subplot
import numpy

class perfect_subplot_group:
    def __init__(self,p_subplot_list,groups=[],get_all=False,logic="or"):
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
        
    def getattrs(self,attr_name,range=False):
        #returns a list of the attribute for all members in the list
        if range:
            # we only return the plotted part of the array attribute
            if attr_name == "data":
                return [p_subplot.data_inrange() for p_subplot in self.p_subplot_list]
            if attr_name == "x":
                return [p_subplot.x_inrange() for p_subplot in self.p_subplot_list]
            if attr_name == "y":
                return [p_subplot.y_inrange() for p_subplot in self.p_subplot_list]
        else:
            # we return the entire array, disregarding xlimits and ylimits
            return [getattr(p_subplot,attr_name) for p_subplot in self.p_subplot_list]


    def get_max(self,attr_name,margin=0):
        #get the maximum values for all quantites in all subplots in this group
        maxes=[]
        for subplot in self.getattrs(attr_name):
            for simulation in subplot:
                maxes=maxes+[numpy.max(simulation)]
        max=numpy.max(maxes)                
        if margin == 0:
            return max
        else:
            return max + margin*self.get_range(attr_name)

    def get_min(self,attr_name,margin=0):
        #get the maximum values for all quantites in all subplots in this group
        mines=[]
        for subplot in self.getattrs(attr_name):
            for simulation in subplot:
                mines=mines+[numpy.min(simulation)]
        min = numpy.min(mines)
        if margin == 0:
            return min
        else:
            return min - margin*self.get_range(attr_name)

    def get_range(self,attr_name):
        ret = numpy.fabs(self.get_max(attr_name)-self.get_min(attr_name))
        if numpy.fabs(ret) <= numpy.fabs(self.get_max(attr_name)*1e-10):
            print "Warning: p_subplot_group: get_range: range is effectively zero. Artificially padding to make plot work."
            ret = self.get_max(attr_name)*0.1
        return ret
        
    def set_middle_attr(self,attr_name,value):
        #sets an attribute of the middle element of a list.
        #This is designed for plots where subplots in the same group are close to each other
        #so that it makes sense to only label one of them.
        l=len(self.p_subplot_list)
        i=int(l/2.0)
        setattr(self.p_subplot_list[i],attr_name,value)

    def set_middle_ylabel(self,ylabel):
        #set the ylabel of a connected group to be in the middle
        #will not be sensible for groups that are spread out
        l=len(self.p_subplot_list)
        i=int(l/2.0)
        setattr(self.p_subplot_list[i],"yaxis_label",ylabel)
        if l%2 == 0:
            #even number of subplots: set ylabel just above int(l/2.0) subplot
            setattr(self.p_subplot_list[i],"yaxis_label_y",1.0)
