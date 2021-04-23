#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 17:49:02 2021

@author: tim
"""
import numpy as np
from comparison_fns import dist, mt_to_l

def valid_events(element,idx):
    '''
    Function for determining whether events are still valid, by comparing the indices

    Parameters
    ----------
    element : Event list element
    idx : index of MT to be compared w/

    Returns
    -------
    None.

    '''
    # print(element[1])
    tuple_ = element[1]
    if type(tuple_) == tuple:
        if tuple_[0] == idx or tuple_[1] ==idx:
            return(False)
        else:
            return(True)
    else:
        if tuple_ == idx:
            # print('False')
            return(False)
        else:
            # print('True')
            return(True)

def tri_dist1(p_c,p1,p2):
    '''
    Distance from centre point p_c to line segment p1-p2 using Heron's formula
    for the area

    Parameters
    ----------
    p_c : Centre point
    p1 : start (or end) of line seg
    p2 : end (or start) of line seg

    Returns
    -------
    Perp distance from centre point to line seg

    '''
    d1 = dist(p_c, p1) #one leg
    d2 = dist(p_c, p2) #other leg
    d_p = dist(p1,p2) #line segment length aka base length
    s = (d1+d2+d_p)/2 #semi-perimeter
    B = s*(s-d1)*(s-d2)*(s-d_p)
    h = 0
    if B < 0:
        h =0
    else:
        A = np.sqrt(s*(s-d1)*(s-d2)*(s-d_p)) #heron's formula
        h = 2*A/d_p #height given by A=1/2 d_p h
    return(h)

def tri_dist(p_c,p1,p2):
    '''
    Distance from centre point p_c to line segment p1-p2 using Heron's formula
    for the area

    Parameters
    ----------
    p_c : Centre point
    p1 : start (or end) of line seg
    p2 : end (or start) of line seg

    Returns
    -------
    Perp distance from centre point to line seg

    '''
    p1x, p1y, p2x, p2y = p1[0], p1[1], p2[0], p2[1] #components of each point
    px, py = p_c[0], p_c[1]
    t = (px-p1x)*(p2x-p1x) + (py-p1y)*(p2y-p1y) #(p-p1)dot(p2-p1)
    t = t/((p2x-p1x)**2 + (p2y-p1y)**2) #|p2-p1|^2
    t = min(max(t,0),1) #point on line
    point = [p1x + t*(p2x-p1x)-px, p1y+t*(p2y-p1y)-py]
    h = point[0]**2 + point[1]**2
    h = np.sqrt(h)
    return(h)

def check_region(idx1,idx2,mt_list,t,r):
    '''
    Checks MT2 (number=idx2) will enter within a radius r of MT1 before MT1 reaches said radius
    Also checks if MT1 and MT2 both shrink (thus don't intersect at all). For use in creating new
    MT index lists

    Parameters
    ----------
    idx1 : Number of MT1
    idx2 : Number of MT2
    mt_list : List of current MTs
    t : current time
    r : radius of region
    Returns
    -------
    True if MT2 will enter the region, False otherwise. Not considering when only one MT shrinks

    '''
    l_idx1, l_idx2 = mt_to_l(mt_list,idx1), mt_to_l(mt_list, idx2) #get list indices
    mt1, mt2 = mt_list[l_idx1], mt_list[l_idx2] #get corresponding MTs
    result = False
    grow1, grow2 = mt1.grow, mt2.grow
    if  grow1==False and grow2==False: #if both are shrinking
        0
    elif grow1 and grow2: #both growing
        dt1 = t-mt1.update_t[-1] #time mt2 grows
        dt2 = t-mt2.update_t[-1] #time mt2 grows
        seg2 = mt2.seg #previous segment points
        for i in  range(len(seg2)):
            p1 = [mt1.seg[-1][0] + dt1*np.cos(mt1.angle[-1]),mt1.seg[-1][1] + dt1*np.sin(mt1.angle[-1])] #current point of second MT
            if i == len(seg2)-1: #dynamic ends need to be compared
                p2 = [mt2.seg[-1][0] + dt2*np.cos(mt2.angle[-1]),mt2.seg[-1][1] + dt2*np.sin(mt2.angle[-1])] #current point of second MT
                d = dist(p2, p1) #distance between them
                if d<2*r: #possibly enters region
                    result = True
            else: #previously established segments
                d = tri_dist(p1, seg2[i], seg2[i+1])
                if d<r: #possibly enters region
                    result = True
            if result: #true stop after we know it's intersected
                break
    elif grow1==False or grow2==False: #one not growing
        if grow1==False:
            mts,mtg = mt1, mt2 #denote which is the shrinking MT
        else:
            mts,mtg = mt2, mt1
        dt1 = mtg.region_t[-1]-mtg.update_t[-1] #time mt2 grows
        seg2 = mts.seg #previous segment points
        for i in  range(len(seg2)-1):
            p1 = [mtg.seg[-1][0] + dt1*np.cos(mtg.angle[-1]),mtg.seg[-1][1] + dt1*np.sin(mtg.angle[-1])] #current point of second MT
            d = tri_dist(p1, seg2[i], seg2[i+1])
            if d<r: #possibly enters region
                result = True
            if result: #true stop after we know it's intersected
                break
    else: #one shrinking, complicated - may edit later
        result = True
    return(result)
