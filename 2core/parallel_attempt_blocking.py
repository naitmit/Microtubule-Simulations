from mpi4py import MPI
from comparison_fns import mt
import numpy as np
import random as rnd
from time import time
from sim_algs_parallel import update_event_list,update, exit_data, undo_update
import sys
# import matplotlib.pyplot as plt
from time import sleep
from parameters import xdomain, ydomain

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
    print("Process has started:",rank)
    sys.stdout.flush()
else:
    print("Process has started:",rank)
    sys.stdout.flush()



rnd.seed(rank)
N = 3#number of MTs
mt_list = []
for i in range(N):
    mt_list.append(mt(i,rank))
    x,y = rnd.uniform(xdomain[0],xdomain[1]), rnd.uniform(ydomain[0],ydomain[1])
    th = rnd.uniform(0,2*np.pi)
    mt_list[i].seg = [[x,y]]
    mt_list[i].angle = [th]
    mt_list[i].region_t.append(0)

t = 0
start = time()
event_list = []
pevent_list = []#past events
del_mt = [] #deleted MTs
cross_data = None
changed_mt = None
policy = None
rewind = False
i=0
crossings = 0
L = 40
conv = L/0.08 #time conversion factor
cross_data = None
dest_ = (rank+1)%2
source_ = dest_
T=600
troubleshoot = 100000 #when to start printing statements for touble shooting

start = MPI.Wtime()
t1 = MPI.Wtime()
t2 = MPI.Wtime()

if rank == 0:
    print('Simulating on', xdomain,'x',ydomain, 'with real length (microm)', L)
    print('Simulating until (s):', T)
    print('Processor grid:', [1,2],'\n')
    sys.stdout.flush()

i1, i1 = 0,0 #for timing
while t*conv < T:
    if rank==0:
        i2 = np.floor(t*conv/10) #prints every ten seconds simulated
        if i2 != i1 and i2>i1: #needed to include the rewinding
            t2 = MPI.Wtime()
            print('Rank 0 at real time (s)', t*conv)
            print('Elapsed wall time (s):', t2-t1)
            print('Length of previous event list: ', len(pevent_list))
            print('Length of event list: ', len(event_list),'\n')
            sys.stdout.flush()
            t1 = MPI.Wtime()
            i1 = i2
    next_event = update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,cross_data,rewind)
    rewind= False
    # if Result[3]*conv > T:
    #     print('Rank', rank, 'plotting')
    #     sys.stdout.flush()
    #     plot_snap(rank, mt_list, t,i,dest='./', save=True)
    update_return = update(rank,N,mt_list,next_event[0],next_event[1],next_event[2], next_event[3], next_event[4],t,cross_data)
    changed_mt = update_return[0]
    policy = next_event[2]
    t = next_event[3]
    i+=1
    '''///Sending Procedure///'''
    if update_return[1] or t>=T: #inform of t>T or crossing events
        if t>=T:
            cross_data = [t]
        else:
            cross_data = exit_data(rank, policy, changed_mt, mt_list,t) #get cross data
        if rank != 0:
            if (t> troubleshoot):
                print('Rank', rank, 'has crossed at time', t,'sending info now')
                sys.stdout.flush()
            # if (t> troubleshoot):
            #     print('Policy before sent',policy)
            #     sys.stdout.flush()

            comm.send(cross_data,dest=0, tag=0)
            # print('Rank', rank, 'waiting for confirmation')
            confirm = comm.recv(source=0,tag=1) #get confirmation that other process can catch up
            if len(confirm) == 1: #only when both processes are done
                print('Rank', rank, 'has recieved confirmation of simluation completion at step ',i,'. Shutting down.')
                sys.stdout.flush()
                break
            if (t> troubleshoot):
                print('Rank', rank, 'has recieved confirmation: ',confirm[0])
                sys.stdout.flush()
            # if (rank == 0) and (t> troubleshoot):
            #     print('Policy before after confirmation',policy)
            #     sys.stdout.flush()
            if not confirm[0]: #if cross happens
                rcross_data = confirm[-1]
                rec_t = rcross_data[1].update_t[-1]
                if (t> troubleshoot):
                    print('Rank', rank, 'will rewind to before', rec_t)
                    sys.stdout.flush()
                while t > rec_t:
                    i-=1
                    t = undo_update(rank, mt_list, pevent_list, event_list, del_mt) #undo events
                rewind = True #clear event list
                update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,rcross_data,rewind) #update list with previous changed MT
                rewind = False #new list generated
                policy = 'comm' #policy will be used to add the new crossing event to list
                cross_data = rcross_data #cross data is from other process
                if (t> troubleshoot):
                    print('Rank', rank, 'has unwound to time:', t)
                    sys.stdout.flush()
            else: #continue as usual
                event_list[:] = [x for x in event_list if x[2] not in ['cross_in','shrink_in']] #there should not be crossings yet
            # if (rank == 0) and (t> troubleshoot):
            #     print('Policy before after confirmation',policy)
            #     sys.stdout.flush()
        elif rank ==0: #root
            if t > troubleshoot:
                print('Rank', rank, 'has crossed at time', t,',waiting to compare')
                sys.stdout.flush()
            rcross_data = comm.recv(source=1, tag=0) #wait to recieve other cross data
            rec_t = None #time of other processor
            if len(rcross_data) == 1: #if other process has passed T
                rec_t = rcross_data[0]
                if t >=T: #both ahead, finish
                    send_confirm = ['done']
                    print('Rank', rank, 'has finished along with rank 1 at step ',i,'. Telling rank 1 to finish.')
                    sys.stdout.flush()
                    comm.send(send_confirm,dest=1,tag=1) #tell other process to finish
                    break
                else: #this process crossed before T
                    if t > troubleshoot:
                        print('Rank', rank, 'crossed first, sending to other rank')
                        sys.stdout.flush()
                    send_confirm = [False,cross_data] #send confirmation data
                    comm.send(send_confirm,dest=1,tag=1)
                    event_list[:] = [x for x in event_list if x[2] not in ['cross_in','shrink_in']] #there should not be crossings yet
            else: #other process is behind T
                rec_t = rcross_data[1].update_t[-1] #other crossing time
                if t<rec_t: #if other process is ahead but both behind T
                    if t > troubleshoot:
                        print('Rank', rank, 'crossed first, sending to other rank')
                        sys.stdout.flush()
                    send_confirm = [False,cross_data] #send confirmation data
                    comm.send(send_confirm,dest=1,tag=1)
                    event_list[:] = [x for x in event_list if x[2] not in ['cross_in','shrink_in']] #there should not be crossings yet
                else:#it's ahead, need to rewind regardless of whether this process is ahead of T
                    if t > troubleshoot:
                        print('Rank', rank, 'crossed is ahead, need to rewind to', rec_t)
                        sys.stdout.flush()

                    while t > rec_t:
                        i-=1
                        t = undo_update(rank, mt_list, pevent_list, event_list, del_mt) #undo events
                    rewind = True #clear event list
                    update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,rcross_data,rewind) #update list with previous changed MT
                    rewind = False #new list generated
                    policy = 'comm' #policy will be used to add the new crossing event to list
                    cross_data = rcross_data #cross_data is from other process
                    if rcross_data[0] in ['new_left','new_right']:
                        assert rcross_data[1].prev_t is not None
                    if t > troubleshoot:
                        print('Rank', rank, 'has unwound to time:', t)
                        sys.stdout.flush()

                    send_confirm = [True,None] #send confirmation data
                    comm.send(send_confirm,dest=1,tag=1)
    if policy in ['nucleate', 'top','bottom','cross_in']:
        N+=1
end = MPI.Wtime()
elapsed_time = end-start
# print(rank, elapsed_time)
# sys.stdout.flush()
total = comm.reduce(elapsed_time, op=max, root=0)
if rank == 0:
    print('Total time taken:', total)
sys.stdout.flush()
