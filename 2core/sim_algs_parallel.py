#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 23:12:47 2021

@author: tim
Simulation algorthims, with parallel edits for 2 processor case
"""
# from mpi4py import MPI
import numpy as np
import itertools as it
from comparison_fns import compare, inter_bdry, dist, which_seg, mt, mt_to_l
from parameters import v_s, L, xdomain, ydomain, tip_scaled
from event_fns import valid_events, check_region
from zippering import zip_cat
import random as rnd
from time import time
from plotting import plot_snap

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
    # k=1
    r, x, y = rnd.uniform(0,1), rnd.uniform(xdomain[0],xdomain[1]), rnd.uniform(ydomain[0],ydomain[1])
    p = [x,y] #point
    dt = np.log(1/r)/k #new time
    T = t+dt
    # T = np.infty
    return(T, p)


def update_event_list(rank, MT_list,event_list,t,pevent_list,del_mt,last_result,policy,cross_data, rewind):
    '''
    Parameters
    ----------
    MT_list : List of MT objects
    pair_list : List of tuples of pairwise MT indices. Each tuple represents MTs to be compared
    t: current time
    pevet_list: list of previous events which occured
    del_mt:
    cross_data: for when MT crosses into processor
    rewind: if true, need to recalculate event_list from existing MTs

    Returns
    -------
    Collision policy
    BE CAREFUL: THE DISTANCE RETURNED IN BDRY COLLISION IS DISTANCE OF PRESENT TIP TO BDRY, NOT
    SEGMENT LENGTH!!!
    '''
    r = tip_scaled*1.1 #radius
    mt_index = [mt.number for mt in MT_list if mt.exist] #indices of active mts
    n = len(mt_index)
    l = None
    # print(t)
    pair_list = None
    if t==0 or rewind: #calculate all possible pairs
        if rewind:
            event_list[:] = [x for x in event_list if x[2]=='nucleate'] #find all nucleation events
        pair_list = list(it.combinations(mt_index,2)) #list of combinations of indices of MTs
        l = len(pair_list)

        for j in range(n): #check for bdry collisions
            MT = MT_list[j]
            if MT.grow == True: #no bdry collision
                bdry_res = inter_bdry(rank,MT, MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,MT.number,bdry_res[2],next_time,bdry_res[1]])
                #regions
                event_list.append([False,MT.number,'new_region',t+r,None])#new when there's a new region
            elif MT.grow ==False and MT.hit_bdry==False: #disappearing events need to be recalculated
                mt_l = np.sum(MT.seg_dist)#total length of shrinking mt, cannot got lower than this
                disap_t = mt_l/v_s+MT.update_t[-1] #collision distances from 1 to 2, growing can only collide w/ shrinking
                event_list.append([False,MT.number,'disap',disap_t,None]) #disappearing event
        if t== 0: #only time we need to new nucleation
            th = rnd.uniform(0,2*np.pi)
            nucleate  = nuc_t(t) #nucleation event
            event_list.append([False,None,'nucleate',nucleate[0],[nucleate[1],th]])
    else:
        if policy == 'comm': #if there's been a communication
            new_policy = cross_data[0] #did an MT shrink or appear?
            mt_cross = cross_data[1] #MT of interest
            if new_policy =='disap': #shrank
                event_list.append([False,None,'shrink_in',mt_cross.disap_t,None]) #shrinking in the future
            else:
                event_list.append([False,None,'cross_in',mt_cross.update_t[-1],None]) #crossing in
        else:
            l_idx = mt_to_l(MT_list, last_result)
            if policy in ['top','bottom']:#if collided w/ bdry, two MTs are updated
                old_idx = MT_list[l_idx].prev_mt #index of previous MT
                new_idx1 = [last_result]
                new_idx2 = [old_idx] #index of updated MTs
                mt_index1= [p for p in mt_index if ((p != last_result) \
                and check_region(last_result,p,MT_list,t,r))] #for comparing new MT
                mt_index2= [p for p in mt_index if ((p != old_idx and p != last_result) \
                or check_region(old_idx,p,MT_list,t,r))] #for comparing old MT
                pair_list1 = list(it.product(mt_index1,new_idx1)) #generate new pairs
                pair_list2 = list(it.product(mt_index2,new_idx2))
                pair_list = pair_list1+pair_list2 #combine tuple lists
                l = len(pair_list)
                event_list[:] = [x for x in event_list if valid_events(x,old_idx)] #discard all invalid events
                #Find bdry events
                MT = MT_list[l_idx]
                bdry_res = inter_bdry(rank, MT, MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,last_result,bdry_res[2],next_time,bdry_res[1]])
                #Find region event
                event_list.append([False,last_result,'new_region',t+r,None])
            elif policy =='cross_in': #crossing in from left or right
                event_list.pop(0)
                new_idx = [last_result]
                mt_index= [p for p in mt_index if ((p != last_result) \
                and check_region(last_result,p,MT_list,t,r))] #for comparing new MT

                pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                l = len(pair_list)
                #Find bdry events
                MT = MT_list[l_idx]
                bdry_res = inter_bdry(rank,MT, MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,last_result,bdry_res[2],next_time,bdry_res[1]])
                #Find region event
                event_list.append([False,last_result,'new_region',t+r,None])
            elif policy in ['left','right']: #these mts have entered a new processor
                '''PARALLEL EDIT'''
                new_idx = [last_result]
                mt_index= [p for p in mt_index if ((p != last_result) \
                and check_region(last_result,p,MT_list,t,r))] #for comparing new MT
                pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                l = len(pair_list)
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                '''PARALLEL EDIT'''
            elif policy == 'deflect': #deflection
                new_idx = [last_result] #index of updated MTs
                mt_index= [p for p in mt_index if p != last_result and check_region(last_result,p,MT_list,t,r)] #indices of new pairs to check over
                n = len(mt_index)
                pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                l = len(pair_list)
                assert(last_result == event_list[0][1])
                assert(valid_events(event_list[0],last_result)==False)
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                MT = MT_list[l_idx]
                #bdry events
                bdry_res = inter_bdry(rank,MT,MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
                #Find region event
                event_list.append([False,last_result,'new_region',t+r,None])
            elif policy == 'disap': #if last one disappeared
                if MT_list[l_idx].from_bdry and MT_list[l_idx].origin == rank: #and it's an extension from the same rank
                    prev_idx = MT_list[l_idx].prev_mt #index of newly shrinking MT
                    event_list[:] = [x for x in event_list if valid_events(x,prev_idx)] #discard newly shrinking one
                    event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                    new_idx = [prev_idx] #new list of indices
                    mt_index= [p for p in mt_index if p != prev_idx and check_region(last_result,p,MT_list,t,r)]
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
            elif policy == 'shrink_in': #shrinking case
                event_list.pop(0)
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                new_idx = [last_result] #new list of indices
                mt_index= [p for p in mt_index if p != last_result and check_region(last_result,p,MT_list,t,r)]
                pair_list = list(it.product(mt_index,new_idx)) #new pairs
                l = len(pair_list)
                #bdry events
                mt_l = np.sum(MT_list[l_idx].seg_dist)#total length of shrinking mt, cannot got lower than this
                disap_t = mt_l/v_s+t #collision distances from 1 to 2, growing can only collide w/ shrinking
                event_list.append([False,last_result,'disap',disap_t,None])
            elif policy in ['1hit2','2hit1']:
                if MT_list[l_idx].grow is False: #if updated MT is due to cat. and it didn't disap
                    new_idx = [last_result] #index of updated MTs
                    mt_index= [p for p in mt_index if p != last_result] #indices of new pairs to check over
                    n = len(mt_index)
                    pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                    l = len(pair_list)
                    event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events

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
                    bdry_res = inter_bdry(rank,MT,MT_list) #find intersection info
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
                MT = MT_list[l_idx]
                bdry_res = inter_bdry(rank,MT, MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
                #new nucleation event
                nuc_list = [x for x in event_list if x[2]=='nucleate']#check current nucleation queue
                if len(nuc_list) == 0: #if there are no mpre, calculate
                    th = rnd.uniform(0,2*np.pi)
                    nucleate  = nuc_t(t) #nucleation event
                    event_list.append([False,None,'nucleate',nucleate[0],[nucleate[1],th]])
                #Find region event
                event_list.append([False,last_result,'new_region',t+r,None])
            elif policy == 'new_region': #need to check new region
                new_idx = [last_result] #index of updated MTs
                mt_index= [p for p in mt_index if p!= last_result and check_region(last_result,p,MT_list,t,r)] #indices of new pairs to check over
                n = len(mt_index)
                pair_list = list(it.product(mt_index,new_idx)) #generate new pairs
                l = len(pair_list)
                event_list[:] = [x for x in event_list if valid_events(x,last_result)] #discard all invalid events
                #bdry events
                MT = MT_list[l_idx]
                bdry_res = inter_bdry(0,MT, MT_list) #find intersection info
                next_time = bdry_res[0] + MT.update_t[-1] #time of collision
                event_list.append([False,new_idx[0],bdry_res[2],next_time,bdry_res[1]])
                #Find region event
                event_list.append([False,last_result,'new_region',t+r,None])
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
    event_list.sort(key=lambda x: x[3]) #sort according to time
    pevent_list.append(event_list[0])#record event occurance
    return(event_list[0])


def update(rank,N,mt_list,collided,idxs,policy,dists,pt,t,cross_data):
    '''
    Updates MT based on next event
    Parameters
    ----------
    rank: rank of process
    N: number of current+past MTs in the system
    mt_list : list of MTs
    collided : bool of whether MTs collided (True) or it was a bdry event
    policy : what type of collision event
    dists : event time
    pt : point of intersection
    idx : tuple of idices if two MTs collide, single index if bdry event
    t : current time
    cross_data: for when an MT crosses into the processor

    Returns
    -------
    Index of updated MT
    Whether event is a crossing

    '''
    crossing = False #whether there has been a crossing
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
        seg_idx = which_seg(mt1,mt2,t) #find segment of intersection
        r = rnd.randint(0, 1)
        zip_res = zip_cat(mt1.angle[-1],mt2.angle[seg_idx],pt,r) #zipper or not
        new_pt,new_angle = zip_res[1], zip_res[0]
        mt1.seg.append(new_pt) #update point
        mt1.update_t.append(dists) #update time
        actual_dist = dist(mt1.seg[-2],mt1.seg[-1])
        mt1.seg_dist.append(actual_dist) #append seg distance
        if not zip_res[2]: #no catastrophe
            mt1.angle.append(new_angle) #TODO: CHANGE ANGLE UPDATE
            mt1.region_t.append(dists) #region change
        else: #catastrophe
            mt1.grow = False
        return(mt_idx1,crossing)
    # if collided is False: #if collision w/ bdry or deflection
    elif policy in ['disap']: #one has disappeared
        mt_idx1 = idxs
        l_idx1 = mt_to_l(mt_list, mt_idx1) #need to convert to indices in the list of existing mts
        assert l_idx1 >= 0
        mt1 =mt_list[l_idx1]
        mt1.exist = False #has disappeared
        mt1.disap_t = dists
        mt1.update_t.append(dists)
        '''PARALLEL EDIT START'''
        if mt1.origin != rank: #it's a continuation of another processor's mt
            assert mt1.from_bdry
            crossing = True #will need to update with crossing_update
            '''PARALLEL EDIT END'''
        elif mt1.from_bdry: #if it's an extension, must update previous MT
            mt_prev_idx = mt1.prev_mt
            l_prev_idx = mt_to_l(mt_list, mt_prev_idx) #convert to list indices
            mt_list[l_prev_idx].grow = False #it must shrink now
            mt_list[l_prev_idx].hit_bdry = False #no longer hitting bdry
            mt_list[l_prev_idx].update_t.append(dists) #update time of new bdry point now shrinking
        return(mt_idx1,crossing)
    elif policy == 'shrink_in':
        mt1 = cross_data[1]
        l_prev_idx  = mt_to_l(mt_list,mt1.prev_mt)
        mt_list[l_prev_idx].grow = False #it must shrink now
        mt_list[l_prev_idx].hit_bdry = False #no longer hitting bdry
        mt_list[l_prev_idx].update_t.append(dists) #update time of new bdry point now shrinking
        return(mt1.prev_mt,crossing)
    elif policy == 'nucleate':
        #TODO: NEED TO THINK ABOUT HOW TO DO THIS BETTER
        mt_idx = N
        mt_list.append(mt(mt_idx,rank)) #introduce new MT
        mt1 = mt_list[-1]
        mt1.update_t[-1] = dists
        mt1.region_t.append(dists) #region change
        mt1.seg = [pt[0]] #assign points and angles
        # th = rnd.uniform(0,2*np.pi)
        mt1.angle = [pt[1]]
        return(mt_idx,crossing)
    elif policy == 'cross_in':
        mt1 = cross_data[1] #corresponding microtuble
        if rank != 0:
            assert mt1.origin != rank
        mt1.number = N #give a number
        mt_list.append(mt1)
        return(N, crossing)
    elif policy =='deflect':
        l_idx = mt_to_l(mt_list,idxs)#convert to list index

        mt1 = mt_list[l_idx] #mt to be updated
        assert l_idx>=0
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
        return(mt1.number,crossing)
    elif policy =='new_region': #new region, no MT change
        l_idx = mt_to_l(mt_list,idxs)
        mt_list[l_idx].region_t.append(dists) #new region calculated at this time
        return(idxs,crossing)
    else: #hit bdry, nucleates new MT
        assert policy in ['top','bottom','left','right']
        l_idx = mt_to_l(mt_list, idxs) #convert to list index
        mt1 = mt_list[l_idx] #mt to be updated
        mt1.seg.append(pt) #update point
        th_og = mt1.angle[-1]
        mt1.angle.append(th_og)
        mt1.update_t.append(dists) #update time
        mt1.grow = False #no longer considered growing
        mt1.hit_bdry = True #hit bdry
        mt_idx = None
        px, py = pt[0], pt[1]
        segment_dist = dist(mt1.seg[-2],mt1.seg[-1])
        mt1.seg_dist.append(segment_dist) #segment length
        assert mt1.checkd()
        if policy in ['top', 'bottom']:
            mt_idx = N #new border mt number
            mt1.ext_mt = mt_idx #let it know which one is its continuation
            mt_list.append(mt(mt_idx,rank)) #add new mt to rep continuation

            mt2 = mt_list[-1]
            mt2.number = mt_idx
            mt2.origin = rank
            mt2.from_bdry = True
            mt2.prev_mt = idxs
            mt2.update_t[-1] = dists #give update time
            mt2.angle.append(th_og) #update angle

            if policy == 'top': #for each case, continue
                assert py == 1
                mt2.seg.append([px,ydomain[0]])
            elif policy == 'bottom':
                assert py==0
                mt2.seg.append([px,ydomain[1]])
            mt2.region_t.append(dists) #region change
            '''PARALLEL EDIT START'''
        else: #crossing!
            assert policy in ['left', 'right']
            mt_idx = idxs #original MT
            crossing = True
            '''PARALLEL EDIT END'''
        assert mt1.checkd()
        return(mt_idx,crossing)

def undo_update(rank, mt_list, pevent_list, event_list, del_mt):
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
    idxs, policy = event[1], event[2] #event details

    # if collided is True: #if dynamic ends collided
    assert policy not in ['shrink_in','comm','cross_in']
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
    # if collided is False: #if collision w/ bdry or deflection
    elif policy =='disap': #one has disappeared
        mt_idx1 = idxs
        l_idx1 = mt_to_l(mt_list, mt_idx1) #need to convert to indices in the list of existing mts
        ldel_idx = mt_to_l(del_mt, mt_idx1) #might not exist, check non-existing mts
        mt1 = None
        if l_idx1 >=0: #if in this list
            mt1 =mt_list[l_idx1]
        else: #if in deleted mts
            assert ldel_idx >= 0
            mt1 = del_mt[ldel_idx]
            del_mt.pop(ldel_idx) #delete from deleted mt list
            mt_list.append(mt1) #re-introduce to mtlist
        assert mt1.number == idxs
        # print(mt1.number)
        mt1.exist = True #has un-dissapeared
        mt1.disap_t = None #undo this
        mt1.update_t.pop(-1)
        if mt1.from_bdry and mt1.origin==rank: #if it's an extension, must un-update previous MT
            mt_prev_idx = mt1.prev_mt
            l_prev_idx = mt_to_l(mt_list, mt_prev_idx) #convert to list indices
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

        event_list.append(event) #reinstate nucleation event
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
        mt1.update_t.pop(-1) #un-update time
        mt1.grow = True #un-no longer considered growing
        mt1.hit_bdry = False #un-hit bdry
        if policy in ['top', 'bottom']:
            mt_idx = mt1.ext_mt #new border mt number
            mt1.ext_mt = None #un-let it know which one is its continuation

            l_idx2 = mt_to_l(mt_list, mt_idx) #idx of continuation
            assert l_idx2 >= 0
            mt_list.pop(l_idx2) #un-add new mt to rep continuation
        mt1.seg_dist.pop(-1) #segment length
    pevent_list.pop(-1) #un-updated the latest event
    t_new = pevent_list[-1][3]
    return(t_new)

def exit_data(rank, policy, idx, mt_list,t):
    '''
    Given that MT w/ number idx in rank is crossing outside, prepare data to be
    sent to other process

    Parameters
    ----------
    t: current time
    rank : rank of process
    policy: left/right or disap
    idx : ID of MT in this rank
    mt_list: list of MTs

    Returns
    -------
    Data

    '''
    data = [] #data to be sent
    if policy == 'disap':
        data.append(policy)
        l_idx = mt_to_l(mt_list, idx) #convert to list index
        mt1 = mt_list[l_idx] #mt that disappered
        data.append(mt1) #index in the other processor
    else:
        assert policy in['right','left']
        if policy=='right': #new event in other process
            data.append('new_left')
        else:
            data.append('new_right')
        l_idx = mt_to_l(mt_list, idx) #convert to list index
        mt1 = mt_list[l_idx] #mt to be updated
        pt = mt1.seg[-1][:] #previous point
        if rank == 0: #if rank 0 (left), find identified edge
            if policy == 'right':
                pt[0] = xdomain[0]
            else:
                pt[0] = xdomain[1]
        else: #rank 1 (right)
            if policy == 'left':
                pt[0] = xdomain[1]
            else:
                pt[0] = xdomain[0]
        th = mt1.angle[-1] #angle
        mt2 = mt(None,rank) #continuation mt
        mt2.angle = [th]
        mt2.seg = [pt]
        mt2.from_bdry = True
        mt2.prev_mt = idx
        mt2.update_t[-1] = t#give update time
        mt2.region_t = [t]
        if len(mt1.seg) > 2: #it changed trajectory before hitting bdry
            mt2.prev_t = mt1.update_t[-2]
        else:
            assert len(mt1.seg)==2
            if mt1.origin == rank: #can search native list for other info
                if mt1.from_bdry:
                    l_idx2 = mt_to_l(mt_list,mt1.prev_mt)#may need to consider ext of ext
                    mtp2 = mt_list[l_idx2]
                    assert len(mtp2.update_t) >= 2 #implicit assumption here that previous mt will be more segs
                        # plot_snap([mt1,mtp2],t,i,'./plots2/',False)
                    mt2.prev_t = mtp2.update_t[-2] #time before it hit wall
                else:
                    mt2.prev_t = mt1.update_t[0]
            else:
                mt2.prev_t = mt1.prev_t #last physical update same as parent MT
        data.append(mt2)
    return(data)


if __name__ == '__main__':
    rank = 0
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
    T = 0.2
    start = time()
    event_list = []
    pevent_list = []#past events
    del_mt = [] #deleted MTs
    cross_data = None
    last_result = None
    policy = None
    rewind = False
    i=0
    crossings = 0
    # up_t = 0
    # mt_t = 0
    conv = L/0.08 #time conversion factor
    cross_data = None
    # while t*conv< 100:
    j = 0
    start = time()
    ti = 0
    tf = 0
    print('Simulating on', xdomain,'x',ydomain)
    while t< T:
        if i%200==0:
            tf = time()
            print('Sim time: ', t)
            print('Wall time elapsed: ',tf-ti)
            # print('Update time: ',up_t)
            # print('MT timL: ', mt_t)
            print('Length of previous event list: ', len(pevent_list))
            print('Length of event list: ', len(event_list),'\n')
            # plot_snap(mt_list,t,i,'./plots2/',False)
            ti = time()
            # up_t = 0
            # mt_t = 0
            pevent_list = []
        # if i >0:
        #     plot_snap(mt_list,t,i,'./plots5/',False)
        # print(t, policy)
        Result = update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, last_result, policy,cross_data,rewind)
        rewind= False
        policy = Result[2]
        update_return = update(rank,N,mt_list,Result[0],Result[1],Result[2], Result[3], Result[4],t,cross_data)
        last_result = update_return[0]
        t = Result[3]
        if update_return[1]: #if there's a crossing
            # print(t)
            # break
            cross_data = exit_data(rank, policy, last_result, mt_list,t) #package data to be transferred
            # i+=1
            # if i > 0:
            #     plot_snap(mt_list,t,i,'./plots5/',False)
            update_event_list(rank, mt_list, event_list,t,pevent_list,del_mt, last_result, policy,cross_data,rewind) #update list with previous changed MT
            policy = 'comm' #policy will be used to add the new event to list
        if policy in ['nucleate', 'top','bottom','cross_in']:
            N+=1
        i+=1
        # if update_return[1]:
        #     while i>1: #undoing stuff
        #         i-= 1
        #         t = undo_update(rank, mt_list, pevent_list)
        #         # if i%100==0:
        #         plot_snap(mt_list,t,i,'./plots5/',False)
        #         rewind=True
        #     j+=1
        # if j ==5:
        #     break
    end = time()
    # print('Percent of crossings: ', crossings/i)
    print('Total time: ',end-start)
