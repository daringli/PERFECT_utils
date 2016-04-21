def sort_species_list(species_list,first=["D","He"],last=["e"]):
    #sorts list so that the first element in "first" that is in species_list will be first, then the second that exists, etc.
    #elements between first and last have no prefered order
    output_list=[]
    if any([f in last for f in first]):
        print "sort_species_list: warning: recieved directive to both place element first and last. List will be extended to place it both first and last."
    
    for f in first:
        if f in species_list:
            output_list.append(f)
    for s in species_list:
        if (s not in last) and (s not in first):
            output_list.append(s)
    #first in last will be last, since list is reversed
    for l in last[::-1]:
        if l in species_list:
            output_list.append(l)

    return output_list

if __name__=="__main__":
    a=["1","2","3","4","5","6"]
    f=["4","3"]
    l=["2","1"]
    print sort_species_list(a,f,l)

    print sort_species_list(["e","D","He","N","C","W"])
    print sort_species_list(["e","He","N"])
