#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 23:12:47 2021

@author: tim
Simulation algorthims
"""
# from mpi4py import MPI
import numpy as np
import itertools as it
from comparison_fns import compare, inter_bdry, dist, which_seg, v_s, mt, mt_to_l
from event_fns import valid_events, check_region
from zippering import zip_cat
from parameters import v_s, L, xdomain, ydomain, tip_scaled
import random as rnd
from time import time
import matplotlib.pyplot as plt
from plotting import plot_snap
rank = 0
def nuc_t(t):
    '''
    Next nucleation time using rate k, from Gillespie

    Parameters
    ----------
    t : Current time
    Returns
    -------
    Next nucleation time and nucleation point

    '''
    # L= 40
    lx = xdomain[1]-xdomain[0]
    ly = ydomain[1]-ydomain[0]
    k = lx*ly*L**3*0.1/8#nucleation rate
    k=1
    r, x, y = rnd.uniform(0,1), rnd.uniform(xdomain[0],xdomain[1]), rnd.uniform(ydomain[0],ydomain[1])
    p = [x,y] #point
    dt = np.log(1/r)/k #new time
    T = t+dt
    # T = np.infty
    return(T, p)

def update_event_list(MT_list,event_list,t,pevent_list,del_mt,last_result = None,policy = None):
    '''
    Parameters
    ----------
    MT_list : List of MT objects
    pair_list : List of tuples of pairwise MT indices. Each tuple represents MTs to be compared
    t: current time
    pevet_list: list of previous events which occured
    del_mt:
    Returns
    -------
    Collision policy
    BE CAREFUL: THE DISTANCE RETURNED IN BDRY COLLISION IS DISTANCE OF PRESENT TIP TO BDRY, NOT
    SEGMENT LENGTH!!!
    '''
    r = 0.15*80/L #radius
    # r=2
    mt_index = [mt.number for mt in MT_list if mt.exist] #indices of active mts
    n = len(mt_index)
    l = None
    # print(t)
    pair_list = None
    if t==0:
        pair_list = list(it.combinations(mt_index,2)) #list of combinations of indices of MTs
        l = len(pair_list)

        for j in range(n): #check for bdry collisions
            if MT_list[j].grow == True: #no bdry collision
                MT = MT_list[j]
                bdry_res = inter_bdry(0,MT, MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,MT.number,bdry_res[2],next_time,bdry_res[1]])
                #regions
                event_list.append([False,MT.number,'new_region',t+r,None])
        nucleate  = nuc_t(t) #nucleation event
        event_list.append([False,None,'nucleate',nucleate[0],nucleate[1]])
    else:
        l_idx = mt_to_l(MT_list, last_result)
        if policy in ['top','bottom', 'left','right']:
        # if MT_list[last_result].from_bdry is True and MT_list[MT_list[last_result].prev_mt].update_t == t: #if collided w/ bdry, two MTs are updated

            old_idx = MT_list[l_idx].prev_mt #index of previous MT
            new_idx1 = [last_result]
            new_idx2 = [old_idx] #index of updated MTs

            mt_index1= [p for p in mt_index if ((p != last_result) \
            and check_region(last_result,p,MT_list,t,r))] #for comparing new MT

            mt_index2= [p for p in mt_index if ((p != old_idx and p != last_result) \
            or check_region(old_idx,p,MT_list,t,r))] #for comparing old MT
            # print(mt_index, last_result)
            # n = len(mt_index)
            # assert n == len(MT_list)-1
            pair_list1 = list(it.product(mt_index1,new_idx1)) #generate new pairs
            pair_list2 = list(it.product(mt_index2,new_idx2))
            pair_list = pair_list1+pair_list2 #combine tuple lists

            l = len(pair_list)
            event_list[:] = [x for x in event_list if valid_events(x,old_idx)] #discard all invalid events
            #Find bdry events
            l_idx = mt_to_l(MT_list, last_result)
            MT = MT_list[l_idx]
            bdry_res = inter_bdry(0,MT, MT_list) #find intersection info
            next_time = bdry_res[0] + MT.update_t[-1] #time of collision
            event_list.append([False,last_result,bdry_res[2],next_time,bdry_res[1]])
            #Find region event
            event_list.append([False,last_result,'new_region',t+r,None])
        elif policy == 'deflect': #deflection
            new_idx = [last_result] #index of updated MTs
            mt_index= [p for p in mt_index if p != last_result and check_region(last_result,p,MT_list,t,r)] #indices of new pairs to check over
            n = len(mt_index)
            pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
            l = len(pair_list)
            event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events

            l_idx = mt_to_l(MT_list, new_idx[0])
            MT = MT_list[l_idx]
            #bdry events
            bdry_res = inter_bdry(0,MT,MT_list) #find intersection info
            next_time = bdry_res[0] + MT.update_t[-1] #time of collision
            event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
            #Find region event
            event_list.append([False,last_result,'new_region',t+r,None])
        elif policy in ['1disap','2disap','disap']:
        # elif MT_list[last_result].exist is False: #if last one disappeared
            # print(2)
            l_idx = mt_to_l(MT_list, last_result)
            if MT_list[l_idx].from_bdry: #and it's an extension
                prev_idx = MT_list[l_idx].prev_mt #index of newly shrinking MT
                event_list[:] = [x for x in event_list if valid_events(x,prev_idx)] #discard newly shrinking one
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events

                new_idx = [prev_idx] #new list of indices
                mt_index= [p for p in mt_index if p != prev_idx]
                pair_list = list(it.product(mt_index,new_idx)) #new pairs
                l = len(pair_list)
                assert MT_list[l_idx].exist==False
                #bdry events
                l_idxp = mt_to_l(MT_list, prev_idx) #find microtuble
                mt_l = np.sum(MT_list[l_idxp].seg_dist)#total length of shrinking mt, cannot got lower than this
                disap_t = mt_l/v_s+t #collision distances from 1 to 2, growing can only collide w/ shrinking
                event_list.append([False,prev_idx,'disap',disap_t,None])

                del_mt.append(MT_list[l_idx]) #add deleted mt to list
                MT_list.pop(l_idx) #remove from existent mt list
            else: #if last one disappeared, no continuation
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                l = 0 #no new comparisons to be made with non-existant MT
                assert MT_list[l_idx].exist==False
                del_mt.append(MT_list[l_idx])
                MT_list.pop(l_idx)
        elif policy in ['1hit2','2hit1']:
            l_idx = mt_to_l(MT_list, last_result)
            if MT_list[l_idx].grow is False: #if updated MT is due to cat. and it didn't disap
                new_idx = [last_result] #index of updated MTs
                mt_index= [p for p in mt_index if p != last_result] #indices of new pairs to check over
                n = len(mt_index)
                pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                l = len(pair_list)
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                # MT = MT_list[new_idx[0]]
                # bdry_res = inter_bdry(MT.seg[-1], MT.angle[-1]) #find intersection info
                # next_time = bdry_res[0] + MT.update_t #time of collision
                # event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])

                mt_l = np.sum(MT_list[l_idx].seg_dist)#total length of shrinking mt, cannot got lower than this
                disap_t = mt_l/v_s+t #collision distances from 1 to 2, growing can only collide w/ shrinking
                event_list.append([False,last_result,'disap',disap_t,None]) #disappearing event
            else: #if two dynamics MTs collide, must be growing
                new_idx = [last_result] #index of updated MTs
                mt_index= [p for p in mt_index if p != last_result and check_region(last_result,p,MT_list,t,r)] #indices of new pairs to check over
                n = len(mt_index)
                pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                l = len(pair_list)
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events

                l_idx = mt_to_l(MT_list, new_idx[0])
                MT = MT_list[l_idx]
                #bdry event
                bdry_res = inter_bdry(0,MT,MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
                #Find region event
                event_list.append([False,last_result,'new_region',t+r,None])
        elif policy == 'nucleate':#nucleation occurs
            event_list.pop(0)
            new_idx = [last_result] #index of updated MTs
            mt_index= [p for p in mt_index if check_region(last_result,p,MT_list,t,r)] #indices of new pairs to check over
            n = len(mt_index)
            # assert n == len(MT_list)-1
            pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
            l = len(pair_list)
            #bdry events
            l_idx = mt_to_l(MT_list, last_result)
            MT = MT_list[l_idx]
            bdry_res = inter_bdry(0,MT, MT_list) #find intersection info
            next_time = bdry_res[0] + MT.update_t[-1] #time of collision
            event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
            #new nucleation event
            nucleate  = nuc_t(t) #nucleation event
            event_list.append([False,None,'nucleate',nucleate[0],nucleate[1]])
            #Find region event
            event_list.append([False,last_result,'new_region',t+r,None])
        elif policy == 'new_region': #need to check new region
            event_list.pop(0)
            new_idx = [last_result] #index of updated MTs
            mt_index= [p for p in mt_index if p!= last_result and check_region(last_result,p,MT_list,t,r)] #indices of new pairs to check over
            n = len(mt_index)
            pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
            l = len(pair_list)
            event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
            #bdry events
            # l_idx = mt_to_l(MT_list, last_result)
            #Find region event
            event_list.append([False,last_result,'new_region',t+r,None])
            
            MT = MT_list[l_idx]
            bdry_res = inter_bdry(0,MT, MT_list) #find intersection info
            next_time = bdry_res[0] + MT.update_t[-1] #time of collision
            event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
    # start = time()
    for i in range(l): #compare all pairs
        mt1_i = pair_list[i][0] #numbers of mts to compare
        mt2_i = pair_list[i][1]
        l_mt = mt_to_l(MT_list,mt1_i,mt2_i) #convert to list index
        l_mt1, l_mt2 = l_mt[0],l_mt[1]
        result = compare(MT_list[l_mt1], MT_list[l_mt2],t) #retrieve output object
        if result.policy not in ['no_collision']: #if there is collision
            assert result.dist >=0 and MT_list[l_mt1].exist and MT_list[l_mt2].exist
            next_time = result.dist+t
            event_list.append([True,pair_list[i],result.policy,next_time,result.point]) #apend event details
    # print(pair_list)
    # print('List time: ', time()-start,l)
    event_list.sort(key=lambda x: x[3]) #sort according to time
    pevent_list.append(event_list[0])#record event occurance
    return(event_list[0])


def update(mt_list,collided,idxs,policy,dists,pt,t):
    '''
    Updates MT based on next event
    Parameters
    ----------
    mt_list : list of MTs
    collided : bool of whether MTs collided (True) or it was a bdry event
    policy : what type of collision event
    dists : event time
    pt : point of intersection
    idx : tuple of idices if two MTs collide, single index if bdry event
    t : current time

    Returns
    -------
    Index of updated MT

    '''
    # pt = np.ndarray.tolist(pt)
    if collided is True: #if dynamic ends collided
        if policy =='1hit2' or policy=='2hit1':
            if policy == '1hit2':
                mt_idx1,mt_idx2 = idxs[0],idxs[1] #index of MT to be updated
            else:
                mt_idx1,mt_idx2 = idxs[1],idxs[0]

            l_idx = mt_to_l(mt_list, mt_idx1,mt_idx2) #convert to list indices
            l_idx1, l_idx2 = l_idx[0],l_idx[1]

            mt1 =mt_list[l_idx1]
            mt2 =mt_list[l_idx2]
            # print(idx1,idx2)
            seg_idx = which_seg(mt1,mt2,t)
            # if mt2.angle[seg_idx] != mt2.angle[-1]:
            #     print(mt2.angle[seg_idx], mt2.angle[-1])
            #     print(mt2.number)
            #     assert mt2.angle[seg_idx] == mt2.angle[-1]
            r = rnd.randint(0, 1)
            zip_res = zip_cat(mt1.angle[-1],mt2.angle[seg_idx],pt,r)# NOT RIGHT, NOT LATEST ANGLE!!
            new_pt,new_angle = zip_res[1], zip_res[0]
            # print(zip_res[0],mt1.angle[-1])
            # new_pt, new_angle = pt, mt1.angle[-1]
            mt1.seg.append(new_pt) #update point
            mt1.update_t.append(dists) #update time
            actual_dist = dist(mt1.seg[-2],mt1.seg[-1])
            mt1.seg_dist.append(actual_dist) #append seg distance
            if not zip_res[2]: #no catastrophe
                mt1.angle.append(new_angle) #TODO: CHANGE ANGLE UPDATE
                mt1.region_t.append(dists) #region change
            else: #catastrophe
                mt1.grow = False
            return(mt_idx1)
    if collided is False: #if collision w/ bdry or deflection
        if policy in ['1disap','2disap','disap']: #one has disappeared
                assert policy =='disap'
                mt_idx1 = idxs
                # elif policy == '1disap':
                #     mt_idx1= idxs[0] #index of MT to be updated
                # else:
                #     mt_idx1= idxs[1]
                # print('disapeared at time', dists,idx1)
                l_idx1 = mt_to_l(mt_list, mt_idx1) #need to convert to indices in the list of existing mts
                mt1 =mt_list[l_idx1]
                # mt1 = mt_list[idx1]
                mt1.exist = False #has dissapeared
                mt1.disap_t = dists
                if mt1.from_bdry: #if it's an extension, must update previous MT
                    mt_prev_idx = mt1.prev_mt
                    l_prev_idx = mt_to_l(mt_list, mt_prev_idx) #convert to list indices
                    # mt_prev = [mt for mt in mt_list if mt.prev_mt==prev_idx][0] #find prev mt
                    mt_list[l_prev_idx].grow = False #it must shrink now
                    mt_list[l_prev_idx].hit_bdry = False #no longer hitting bdry
                    mt_list[l_prev_idx].update_t.append(dists) #update time of new bdry point now shrinking

                    # print('previous index', prev_idx)
                return(mt_idx1)
        elif policy == 'nucleate':
            # print('AA')
            #TODO: NEED TO THINK ABOUT HOW TO DO THIS BETTER
            mt_idx = None
            if len(mt_list)==0 and len(del_mt)==0:
                mt_idx=0
            else:
                mt_idx = mt_list[-1].number + 1 #new mt number
            mt_list.append(mt(mt_idx,0)) #introduce new MT
            mt1 = mt_list[-1]
            mt1.update_t[-1] = dists
            mt1.region_t.append(dists) #region change
            mt1.seg = [pt] #assign points and angles
            th = rnd.uniform(0,2*np.pi)
            mt1.angle = [th]
            return(mt_idx)
        elif policy =='deflect':
            l_idx = mt_to_l(mt_list,idxs)#convert to list index

            mt1 = mt_list[l_idx] #mt to be updated
            mt1.seg.append(pt) #update point
            power = rnd.randint(0,1) #to be used to determine deflection +/-
            dth  =np.pi/16 #change in angle
            th_new = mt1.angle[-1] + (-1)**power*dth #deflection
            if th_new> 2*np.pi:
                th_new -= 2*np.pi
            elif th_new <0:
                th_new += 2*np.pi
            mt1.angle.append(th_new)

            segment_dist = dist(mt1.seg[-2],mt1.seg[-1])
            mt1.seg_dist.append(segment_dist) #segment length

            mt1.update_t.append(dists) #update time
            mt1.region_t.append(dists) #region change
            # print('forward', mt1.number, len(mt1.seg_dist))
            return(mt1.number)
        elif policy =='new_region': #new region, no MT change
            l_idx = mt_to_l(mt_list,idxs)
            mt_list[l_idx].region_t.append(dists) #new region calculated at this time
            return(idxs)
        else: #hit bdry, nucleates new MT
            l_idx = mt_to_l(mt_list, idxs) #convert to list index
            mt1 = mt_list[l_idx] #mt to be updated
            mt1.seg.append(pt) #update point
            th_og = mt1.angle[-1]
            mt1.angle.append(th_og)
            # mt1.prev_t = mt1.update_t[-1] #record previous update time
            mt1.update_t.append(dists) #update time
            mt1.grow = False #no longer considered growing
            mt1.hit_bdry = True #hit bdry

            mt_idx = mt_list[-1].number+1 #new border mt number
            mt1.ext_mt = mt_idx #let it know which one is its continuation
            mt_list.append(mt(mt_idx,0)) #add new mt to rep continuation

            mt2 = mt_list[-1]
            mt2.number = mt_idx
            mt2.from_bdry = True
            mt2.prev_mt = idxs
            mt2.update_t[-1] = dists #give update time
            mt2.angle.append(th_og) #update angle
            px, py = pt[0], pt[1]

            if policy == 'top': #for each case, continue
                # assert py == 1
                mt2.seg.append([px,ydomain[0]])
            elif policy == 'bottom':
                # assert py==0
                mt2.seg.append([px,ydomain[1]])
            elif policy =='left':
                # assert px==0
                mt2.seg.append([xdomain[1],py])
            elif policy=='right':
                # assert px==1
                mt2.seg.append([xdomain[0],py])
            segment_dist = dist(mt1.seg[-2],mt1.seg[-1])
            mt1.seg_dist.append(segment_dist) #segment length
            mt2.region_t.append(dists) #region change
            return(mt_idx)

def undo_update(mt_list, pevent_list):
    '''
    Undoes updates
    Parameters
    ----------
    mt_list : list of MTs
    collided : bool of whether MTs collided (True) or it was a bdry event
    policy : what type of collision event
    dists : event time
    pt : point of intersection
    idx : tuple of idices if two MTs collide, single index if bdry event
    t : current time

    Returns
    -------
    Time of rewinded step

    '''
    #TODO UNDO REGION CHANGE: POP IN EACH ONE MT IS STILL GROWING AND ADD ITS OWN CATEGORY!!
    event = pevent_list[-1] #previous event
    collided, idxs, policy = event[0], event[1], event[2] #event details

    if collided is True: #if dynamic ends collided
        if policy =='1hit2' or policy=='2hit1':
            if policy == '1hit2':
                mt_idx1,mt_idx2 = idxs[0],idxs[1] #index of MT to be updated
            else:
                mt_idx1,mt_idx2 = idxs[1],idxs[0]
            l_idx = mt_to_l(mt_list, mt_idx1,mt_idx2) #convert to list indices
            l_idx1, l_idx2 = l_idx[0],l_idx[1]
            assert l_idx1>=0 and l_idx2>=0
            mt1 =mt_list[l_idx1]

            mt1.seg.pop(-1) #un-update point
            mt1.update_t.pop(-1) #un-update time
            mt1.seg_dist.pop(-1) #append seg distance
            if mt1.grow: #no catastrophe
                mt1.angle.pop(-1) #undo angle change
                mt1.region_t.pop(-1)
            else: #catastrophe
                mt1.grow = True #undo cat
    if collided is False: #if collision w/ bdry or deflection
        if policy in ['1disap','2disap','disap']: #one has disappeared
            if policy =='disap':
                mt_idx1 = idxs
            elif policy == '1disap':
                mt_idx1= idxs[0] #index of MT to be updated
            else:
                mt_idx1= idxs[1]
            l_idx1 = mt_to_l(mt_list, mt_idx1) #need to convert to indices in the list of existing mts
            ldel_idx = mt_to_l(del_mt, mt_idx1) #might not exist, check non-existing mts
            mt1 = None
            if l_idx1 >=0: #if in this list
                mt1 =mt_list[l_idx1]
            else: #if in deleted mts
                assert ldel_idx >= 0
                mt1 = del_mt[ldel_idx]
                mt_list.append(mt1) #re-introduce to mtlist

            mt1.exist = True #has un-dissapeared
            mt1.disap_t = None #undo this
            if mt1.from_bdry: #if it's an extension, must un-update previous MT
                mt_prev_idx = mt1.prev_mt
                l_prev_idx = mt_to_l(mt_list, mt_prev_idx) #convert to list indices
                # mt_prev = [mt for mt in mt_list if mt.prev_mt==prev_idx][0] #find prev mt
                mt_list[l_prev_idx].grow = False #it must shrink now
                mt_list[l_prev_idx].hit_bdry = True #hitting bdry
                mt_list[l_prev_idx].update_t.pop(-1) #un-update time of new bdry point
        elif policy == 'nucleate':
            t = event[3] #in a nucleation event, we do not know the MT index ahead of time
            mt_list2 = [mt for mt in mt_list if mt.update_t[-1]==t] #must find the MT this way
            assert len(mt_list2)==1

            mt_idx = mt_list2[0].number #new mt number
            l_idx = mt_to_l(mt_list, mt_idx) #index in list
            assert l_idx >= 0
            mt_list.pop(l_idx) #un-introduce new MT
        elif policy =='new_region': #new region, no MT change
            l_idx = mt_to_l(mt_list,idxs)
            mt_list[l_idx].region_t.pop(-1) #new region calculated at this time
        elif policy =='deflect':
            l_idx = mt_to_l(mt_list,idxs)#convert to list index
            mt1 = mt_list[l_idx] #mt to be updated
            # print('backward',mt1.number, len(mt1.seg_dist))
            mt1.seg.pop(-1) #un-update point
            mt1.angle.pop(-1)
            mt1.seg_dist.pop(-1) #segment length
            mt1.update_t.pop(-1) #un-update time
            mt1.region_t.pop(-1)
        else: #hit bdry, nucleates new MT
            l_idx = mt_to_l(mt_list, idxs) #convert to list index
            mt1 = mt_list[l_idx] #mt to be updated
            mt1.seg.pop(-1) #update point

            mt1.angle.pop(-1)
            # mt1.prev_t = None #un-record previous update time
            mt1.update_t.pop(-1) #un-update time
            mt1.grow = True #un-no longer considered growing
            mt1.hit_bdry = False #un-hit bdry

            mt_idx = mt1.ext_mt #new border mt number
            mt1.ext_mt = None #un-let it know which one is its continuation

            l_idx2 = mt_to_l(mt_list, mt_idx) #idx of continuation
            assert l_idx2 >= 0
            mt_list.pop(l_idx2) #un-add new mt to rep continuation
            mt1.seg_dist.pop(-1) #segment length
    pevent_list.pop(-1) #un-updated the latest event
    t_new = pevent_list[-1][3]
    return(t_new)


# def plot_snap(mt_list,t,k,dest='./plots3/',save=True):
#     l = len(mt_list)
#     for i in range(l):
#         coord = mt_list[i].seg #get coordinates of segments
#         x_seg = [p[0] for p in coord]
#         y_seg = [p[1] for p in coord]
#         th = mt_list[i].angle[-1] #angle
#         t_i = mt_list[i].update_t[-1] #last update time
#         if mt_list[i].exist: #if it exists, plot segments
#             plt.scatter(x_seg,y_seg,s=2, color='black') #scatter of all segment points
#         if t_i == t: #if it's been updated at this time
#             if mt_list[i].exist:
#                 if mt_list[i].grow or mt_list[i].hit_bdry: #growing or stationary
#                     x_i = x_seg[:-1] #all but most recent update
#                     y_i = y_seg[:-1]
#                     x_f = x_seg[-2:] #updated segment
#                     y_f = y_seg[-2:]
#                     plt.plot(x_i,y_i,color='green',linewidth=0.5)
#                     if len(x_f)>1:
#                         plt.plot(x_f,y_f,color='red',linewidth=3, label='MT event')
#                 else:
#                     plt.plot(x_seg,y_seg,color='green',linewidth=0.5) #started to shrink
#             else:
#                 plt.scatter([x_seg[0]],[y_seg[0]], s=10,color='purple') #yellow points are dissapeared MTs
#         else: #if it's continuing to grow at this point in time
#             if mt_list[i].exist:
#                 if mt_list[i].grow or mt_list[i].hit_bdry: #TODO actually grows beyond
#                     x2 = x_seg[-1] + np.cos(th)*(t-t_i)
#                     y2 = y_seg[-1] + np.sin(th)*(t-t_i)
#                     if x2 > 1: #for plotting the tip, so that it's within [0,1]^2
#                         x2 = 1
#                     elif x2<0:
#                         x2 = 0

#                     if y2 > 1:
#                         y2 = 1
#                     elif y2<0:
#                         y2 = 0

#                     x_f = [x_seg[-1], x2] #where it has grown to
#                     y_f = [y_seg[-1], y2]
#                     plt.plot(x_seg,y_seg,color='green',linewidth=0.5)
#                     plt.plot(x_f,y_f, color='blue', linewidth=0.5)
#                 else: #TODO should be shrinking
#                     # print('plotting MT', i)
#                     assert not mt_list[i].grow
#                     cumsum = np.append([0],np.cumsum(mt_list[i].seg_dist)) #cumulative sums
#                     d = (t- t_i)*v_s #distance traversed
#                     # print('distance shrank:', d)
#                     left_over = cumsum[-1] - d
#                     # print('left over', cumsum[-1])
#                     j = 0
#                     while left_over > cumsum[j]: #find index within cumsum
#                         j+=1
#                         if j > len(cumsum)-1:
#                             break
#                     j -= 1
#                     delta = left_over - cumsum[j] #component sticking out from segment
#                     # print(delta, left_over,cumsum[i])
#                     x_i = x_seg[:j+1] #untouched points in shrinkage
#                     y_i = y_seg[:j+1]
#                     # print(x_i,y_i,i)
#                     # print('mt', i)
#                     # print('seg index', j)
#                     th = mt_list[i].angle[j]
#                     x_f = [x_seg[j], x_seg[j] + np.cos(th)*delta] #where it has grown to
#                     y_f = [y_seg[j], y_seg[j] + np.sin(th)*delta]
#                     plt.plot(x_i,y_i,color='green',linewidth=0.5)
#                     plt.plot(x_f,y_f,color='cyan',linewidth=.5)
#             else:
#                 plt.scatter([x_seg[0]],[y_seg[0]],s=10,color='purple')
#                 # print('sdsd')

    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.title('Time {}, frame {}'.format(t,k))
    # plt.legend()
    plt.gca().set_aspect('equal',adjustable='box')
    figure = plt.gcf()
    figure.set_size_inches(10, 10)
    name = 'plot_'+str(k)
    if save:
        plt.savefig(dest+name)
    plt.show()
    plt.clf()
    return(0)



if __name__ == '__main__':
    rnd.seed(0)
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
    last_result = None
    policy = None
    ti = 0
    tf = 0
    i=0
    crossings = 0
    # up_t = 0
    # mt_t = 0
    conv = L/0.08 #time conversion factor
    # while t*conv< 100:
    for i in range(10):
        # if i%100==0:
        #     # tf = time()
        #     # print('Sim time: ', t*conv)
        #     # print('Wall time elapsed: ',tf-ti)
        #     # # print('Update time: ',up_t)
        #     # # print('MT timL: ', mt_t)
        #     # print('Length of previous event list: ', len(pevent_list))
        #     # print('Length of event list: ', len(event_list),'\n')
        #     # plot_snap(mt_list,t,i,'./plots2/',False)
        #     # ti = time()
        #     # # up_t = 0
        #     # # mt_t = 0
        #     pevent_list = []

        # start_up = time()
        if i >0:
            plot_snap(mt_list,t,i,'./plots5/',False)
        Result = update_event_list(mt_list, event_list,t,pevent_list,del_mt, last_result,policy)
        # end_up = time()
        policy = Result[2]
        if policy in ['left','right']: #how many crossings
            # print('Crossing ocurred')
            crossings+=1
        # elif policy == 'new_region':
        #     print('Changed region',t)
        # start_mt = time()
        last_result = update(mt_list,Result[0],Result[1],Result[2], Result[3], Result[4],t)
        # end_mt = time()
        # print(last_result)
        # mt_t += end_mt-start_mt
        # up_t += end_up - start_up
        t = Result[3]
        # print(t)
        i+=1
    # while t>0:
    #     i-= 1
    #     t = undo_update(mt_list, pevent_list)
    #     if i%100==0:
    #         plot_snap(mt_list,t,i,'./plots5/',False)

    end = time()
    print('Percent of crossings: ', crossings/i)
    print('Total time: ',end-start)
