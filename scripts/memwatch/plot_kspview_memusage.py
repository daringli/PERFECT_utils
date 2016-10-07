#!/usr/bin/python

import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm



def data_from_output(output_filename):
    #get memory usage from perfect outputs with -ksp_view
    file = open(output_filename, 'rb')

    record = False
    memory = []
    proc = []

    #read data from output
    for line in file:
        #Data follows line matching this
        if line.strip() == "INFO(16) (size of (in MB) MUMPS internal data used during numerical factorization):":
            record = True
            continue
        
        if record == True:
            l = line.strip()
            if l[0]=="[":
                temp = l.split('[')[1].rsplit(']')
                proc.append(int(temp[0]))
                memory.append(int(temp[1]))
            else:
                #No more data present
                record = False
    #typically, an output can contain several simulations.
    #we split output from several simulations by looking for the zero processor
    proc_lists = [] #list of list of procs in simulations
    mem_lists = [] #list of list of proc memory usage in simulations

    this_proc_list = []
    this_mem_list = []
    p_old = -1
    for p,m in zip(proc,memory):
        if p <= p_old:
            #we have jumped back to zero
            proc_lists.append(numpy.array(this_proc_list))
            mem_lists.append(numpy.array(this_mem_list))
            this_proc_list = []
            this_mem_list = []
        this_proc_list.append(p)
        this_mem_list.append(m)
        p_old = p
    #add the final proc_list
    proc_lists.append(numpy.array(this_proc_list))
    mem_lists.append(numpy.array(this_mem_list))

    # make a bar plot out of memory usage
    N_simuls = len(proc_lists)
    fig, ax = plt.subplots()
    width = 1.0/(N_simuls + 1)
    colorlist = cm.rainbow(numpy.linspace(0,1,N_simuls))
    for i in range(N_simuls):
        x=numpy.arange(len(mem_lists[i]))
        avg_mem = sum(mem_lists[i])/len(mem_lists[i])
        sd_mem = numpy.sqrt(sum((mem_lists[i] - avg_mem)**2)/len(mem_lists[i]))
        rects = ax.bar(x + i*width, mem_lists[i], width, color=colorlist[i])
        ax.axhline(y = avg_mem, color = colorlist[i], linestyle='solid')
        ax.set_xticks(x + (i+1)*width/2.0)
        ax.set_xticklabels(proc_lists[i])
        #autolabel(ax,rects)
        print "simulation: "+ str(i)
        print "avg memory usage: "+ str(avg_mem)
        print "standard deviation (not of mean): "+ str(sd_mem)
        print "sd/avg: "+ str(sd_mem/avg_mem)
        
    plt.show()


def autolabel(ax,rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')
        
if __name__ == "__main__":
    output_filename=sys.argv[1]
    data_from_output(output_filename)
        
        
        
            
