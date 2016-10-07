#!/usr/bin/python

import subprocess
import sys
import time
import argparse


def memwatch(pid_list,help=False,dt=False):
    """monitor virtual memory usage of processes given by the PID list
    if help=True, prints header explaining different columns
    if dt=number, prints output every dt seconds"""
    

    #pids should be strings
    for pid in pid_list:
        pid=str(pid)

    N = len(pid_list)
    
    vm_command="vmstat"
    process_m_command="pmap"


    if help:
        if dt == False:
            print "#Memory usage (KB)"
        else:
            print "#Memory usage (KB), every " + str(dt) + " seconds"
        pid_list_string="".join(["PID:"+pid+' , ' if i is not N-1 else "PID:"+pid for i,pid in enumerate(pid_list)])
        print "#Total VM  , " + pid_list_string
    while True:
        virtual_memory_used=subprocess.check_output([vm_command]).split('\n')[2].split()[2]

        process_memory_usage_list=[]
        for pid in pid_list:
            if process_m_command == "pmap":
                process_memory_usage_list.append(subprocess.check_output([process_m_command,pid]).split('\n')[-2].split()[1][0:-1])
            else:
                print "MEMWATCH ERROR: Unsupported process memory command. Exiting..."
                sys.exit(1)
    
        print str(virtual_memory_used) + ", " + "".join([str(pmu)+', ' if i is not N-1 else str(pmu) for i,pmu in enumerate(process_memory_usage_list)])
        if dt == False:
            break;
        else:
            time.sleep(dt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Monitor total virtual memory usage and memory usage of process(es) given by one or more PIDs',
                                     epilog='Example: memwatch -t <delay between watches in seconds> PID1 PID2 ...')
    parser.add_argument('PIDs', metavar='PID', type=str, nargs='+',
                    help='PID of process to monitor')
    parser.add_argument("-t","--dt",type=int,help='Time between updates. If unspecified, only print once',default=False)
    parser.add_argument("--nohead",help='Do not display table header', action="store_true")

    args = parser.parse_args()
    pid_list= args.PIDs
    try:
        memwatch(pid_list,help=(args.nohead == False),dt=args.dt)
    except KeyboardInterrupt:
        sys.exit(0)

