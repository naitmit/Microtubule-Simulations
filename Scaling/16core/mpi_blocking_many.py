from mpi4py import MPI
from comparison_fns import mt
import numpy as np
import random as rnd
from time import time
from sim_algs_parallel_many import update_event_list,update, exit_data, undo_update
from plotting import plot_snap
import sys
from parameters import xdomain, ydomain, grid, L

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

pevent_list = []#past events
def decide_thread(grid, data, T, size):
    '''
    Given the data from all (>2) processors crossing, decide the earliest
    '''
    return_data = None
    l, w = grid[0], grid[1] #length, width of threads
    times = [] #crossing times of processors
    for i in range(l*w): #iterate through data and record times
        data_i = data[i]
        if len(data_i) == 1: #processor has reached end time
            t_i = data_i[0]
        else:
            t_i = data_i[1].update_t[-1]
        times.append(t_i) #record time
    r = np.argmin(times) #find minimum time processor
    if times[r] >= T: #all processors have completed
        return_data = [times[r]]
    else:
        r_rec = None #rank to be communicated with
        t_cross = times[r]
        data_r = data[r]
        cross_policy = data_r[0] #find policy
        m = np.floor(r/w) #which row the process lies on
        if cross_policy == 'disap': #shrinking in
            r_rec = data_r[1].origin
        elif cross_policy == 'new_left': #rank r hit the right
            r_rec = (r-m*w +1)%w + m*w
        elif cross_policy == 'new_right': #rank r hit the left
            r_rec = (r-m*w -1)%w + m*w
        elif cross_policy == 'new_top': #rank r hit the bottom
            r_rec = (r-w)%size
        elif cross_policy == 'new_bottom': #rank r hit the top
            r_rec = (r+w)%size
        return_data = [[r,r_rec], data_r] #data to be returned
        # print(cross_policy)
        if l>1:
            assert r != r_rec
        # print('Recieving rank', r_rec)
    return(return_data)

print("Process has started:",rank)
sys.stdout.flush()

rnd.seed(rank)
N = 1#number of MTs
mt_list = []
for i in range(N):
    mt_list.append(mt(i,rank))
    x,y = rnd.uniform(xdomain[0],xdomain[1]), rnd.uniform(ydomain[0],ydomain[1])
    th = rnd.uniform(0,2*np.pi)
    mt_list[i].seg = [[x,y]]
    mt_list[i].angle = [th]
    mt_list[i].region_t.append(0)
    # pevent_list. append([False,None,'nucleate',0,[[x,y],th]])

t = 0
start = time()
event_list = []

del_mt = [] #deleted MTs
cross_data = None
changed_mt = None
policy = None
next_event = None
rewind = False
i=0
crossings = 0
conv = L/0.08 #convert simulation time to seconds
cross_data = None

T=600
troubleshoot = T+1 #when to start printing statements for touble shooting

start = MPI.Wtime()
t1 = MPI.Wtime() #timing between iteration chunks
t2 = MPI.Wtime()
i1, i1 = 0,0 #for timing

cross = 0

if rank == 0:
    print('Simulating on', xdomain,'x',ydomain, 'with real length (microm)', L)
    print('Simulating until (s):', T)
    print('Processor grid:', grid,'\n')
    sys.stdout.flush()

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
        # plot_snap(mt_list,t,i,'./rank'+str(rank)+'/')
    next_event = update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,cross_data,rewind)
    rewind= False
    if next_event[3]*conv > T:
        print('Rank', rank, 'plotting')
        sys.stdout.flush()
        plot_snap(rank, mt_list, t,i,dest='./', save=True)
    update_return = update(rank,N,mt_list,next_event[0],next_event[1],next_event[2], next_event[3], next_event[4],t,cross_data)
    changed_mt = update_return[0]
    policy = next_event[2]
    t = next_event[3]
    i+=1
    '''///Sending Procedure///'''
    if update_return[1] or t>=T: #inform of t>T or crossing events
        cross+=1 #add crossing
        if t>=T:
            cross_data = [t]
        else:
            cross_data = exit_data(rank, policy, changed_mt, mt_list,t) #get cross data
        if (t> troubleshoot):
            print('Rank', rank, 'has crossed at time', t,'sending info now')
            sys.stdout.flush()
        '''PARALLEL EDIT'''
        data = None
        data = comm.gather(cross_data, root=0) #gather crossover info into root
        rcross_data = None
        if rank == 0: #comaprison happens on root process
            rcross_data = decide_thread(grid, data, T, size) #find the earliest and relevant information
        rcross_data = comm.bcast(rcross_data, root=0) #broadcast the results
        '''             '''
        rec_t = None #recieved time
        if len(rcross_data) == 1: #only when both processes are done
            rec_t = rcross_data[0]
            print('Rank', rank, 'has recieved confirmation of simluation completion at step ',i,' at time ', t)
            sys.stdout.flush()
            break
        else: #N-1 processors need to rewind
            rec_t = rcross_data[1][1].update_t[-1] #crossing time
            r, r_rec = rcross_data[0][0], rcross_data[0][1]#rank of where the crossing event was calculated and rank being crossed into
            if rank == r_rec and rank == r: #if it crosses into itself
                if (t> troubleshoot):
                    print('Rank', rank, 'had MT cross back into it')
                    sys.stdout.flush()
                update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,rcross_data,rewind) #update list with previous changed MT
                rewind = False #new list generated
                policy = 'comm' #policy will be used to add the new crossing event to list
                cross_data = rcross_data[1] #cross data is from other process
                if (t> troubleshoot):
                    print('Crossing back in with data', cross_data)
                    sys.stdout.flush()
            elif rank == r: #continue as usual
                if t > troubleshoot:
                    print('Rank', rank, 'crossed first, proceeding as usual')
                    sys.stdout.flush()
                event_list[:] = [x for x in event_list if x[2] not in ['cross_in','shrink_in']] #there should not be crossings yet
                # update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,cross_data,rewind)
            elif rank == r_rec:
                if (t> troubleshoot):
                    print('Rank', rank, 'will rewind to before', rec_t)
                    sys.stdout.flush()
                while t > rec_t:
                    i-=1
                    t = undo_update(rank, mt_list, pevent_list, event_list, del_mt) #undo events
                if (t> troubleshoot):
                    print('Rank', rank, 'has unwound to time:', t)
                    sys.stdout.flush()
                rewind = True #clear event list
                if (t> troubleshoot):
                    print('Rank', rank, 'had MT cross into it')
                    sys.stdout.flush()
                update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, changed_mt, policy,rcross_data,rewind) #update list with previous changed MT
                rewind = False #new list generated
                policy = 'comm' #policy will be used to add the new crossing event to list
                cross_data = rcross_data[1] #cross data is from other process
            else:
                 if (t> troubleshoot):
                    print('Rank', rank, 'will rewind to before', rec_t)
                    sys.stdout.flush()
                 while t > rec_t:
                    i-=1
                    t = undo_update(rank, mt_list, pevent_list, event_list, del_mt) #undo events
                    sys.stdout.flush()
                 if (t> troubleshoot):
                    print('Rank', rank, 'has unwound to time:', t)
                    sys.stdout.flush()
                 rewind = True #clear event list
            # if (rank == 0) and (t> troubleshoot):
            #     print('Policy before after confirmation',policy)
            #     sys.stdout.flush()
    if policy in ['nucleate', 'top','bottom','cross_in']:
        N+=1
end = MPI.Wtime()
elapsed_time = end-start
total = comm.reduce(elapsed_time, op=max, root=0)
if rank == 0:
    print('Number of events:', i)
    print('Number of crossings:', j)
    print('Total time taken:', total)
    sys.stdout.flush()
