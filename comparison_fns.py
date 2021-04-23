#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 15:14:24 2021

@author: tim

For functions used in comparing two points
"""
import numpy as np
from parameters import v_s, xdomain, ydomain, tip_scaled

class mt:
    """
    Class for indiv. MTs
    """
    def __init__(self,number,origin):
        self.number = number
        self.exist = True #whether MT has disapeared due to shrinking
        self.seg = [] #nx2 array to store previous segment points
        self.seg_dist = [] #n-1 array to store distance. ith is dist. from seg[i] to seg[i+1]
        self.angle = []# n array angles of trajectories correpsonding to seg points

        self.update_t = [0.0] #time of last segment point update

        self.grow = True #bool for growth or shrink

        self.hit_bdry = False #bool for whether it hit bdry
        self.ext_mt = None #int for MT index of continued MT

        self.from_bdry = False #bool for whether it is extension from an MT which hit bdry
        self.prev_mt = None #MT index from which it is an exension of
        self.disap_t = None #time of disappearance
        self.prev_t = None #previous trajectory change time if it's a continuation
        self.tip_l = tip_scaled #tip length
        # self.rank = None #rank from which this originated
        self.region_t = [] #when the last region calculated
        self.origin = origin #original ID: [rank, ID in that rank]
    def checkd(self):
        tf = True
        for i in range(len(self.seg_dist)):
            if dist(self.seg[i],self.seg[i+1]) !=self.seg_dist[i]:
                tf = False
        return(tf)

class compare_return:
    '''
    Class for return values of compare function
    '''
    def __init__(self):
        self.policy = None #string: no_collision, 1hit2, 2hit1, 1disap, 2disap, hit_bdry
        self.dist = None #float for collision distance
        self.point = None #point of collision

def mt_to_l(mt_list,mt_idx1,mt_idx2=-1):
    l_idx1, l_idx2 = -1,-1 #need to convert to indices in the list of existing mts
    if mt_idx2>=0:
        asgn1,asgn2 = False,False #bools for whether assigned
        for i in range(len(mt_list)):
            if mt_list[i].number == mt_idx1:
                l_idx1 = i
                asgn1 = True
            if mt_list[i].number== mt_idx2:
                l_idx2 = i
                asgn2 = True
            if asgn2 and asgn1:
                break #no more looping
        return(l_idx1,l_idx2)
    else:
        for i in range(len(mt_list)):
            if mt_list[i].number == mt_idx1:
                l_idx1 = i
                break
        return(l_idx1)

def dist(p1,p2,square=False):
    '''
    Parameters
    ----------
    p1 : 2D array
    p2 : 2D array
    square : Bool, optional
        Decide whether to return sqaure distance. The default is False.

    Returns
    -------
    Distance between p1,p2 as float
    '''
    dx = p1[0] - p2[0] #calculate distances
    dy = p1[1] - p2[1]
    d = dx**2 + dy**2
    if square is False:
        return(np.sqrt(d))
    else:
        return(d)

def inter(p1,p2,theta1,theta2):
    '''
    Parameters
    ----------
    p1 : 2D array
    p2 : 2D array
    theta1 : Float of first angle, slope
    theta2 : Float of second angle, slope

    Returns
    -------
    Bool for intersection
    2D array of intersection between two lines
    2D array of distances from p1, p2 respectively
    '''
    assert theta1 != np.pi/2 and theta2 != np.pi/2 and theta1 != 3*np.pi/2 and theta2 != 3*np.pi/2

    if theta1 == theta2:
        return(False, None, None)
    else:
        m1,m2 = np.tan(theta1), np.tan(theta2) #slope intercept form
        if m1==m2:
            return(False, None, None)
        else:
            x1,x2,y1,y2 = p1[0],p2[0],p1[1],p2[1]
            b1,b2 = y1-x1*m1, y2-x2*m2
            x = (b2-b1)/(m1-m2) #equations for intersections
            y = m1*(b2-b1)/(m1-m2) + b1

            P = [x,y] #points and distances
            D = [dist(P,p1),dist(P,p2)]
            # print(P)
            if np.sign(x-x1)==np.sign(np.cos(theta1)) and np.sign(x-x2)==np.sign(np.cos(theta2)): #if in trajectory
                return(True, P, D)
            else:
                return(False, None, None)

def inter_bdry(rank,MT, MT_list):
    '''
    Parameters
    ----------
    MT: MT class
    Returns
    -------
    Float. Distance intersection with boundry
    Array. Point of intersection
    String. Which wall it hits
    '''
    p, th = MT.seg[-1], MT.angle[-1]
    xi, xf = xdomain[0],xdomain[1] #start, end of x domain
    yi, yf = ydomain[0],ydomain[1] #start, end of y domain
    x,y = p[0],p[1] #coordinates of points

    P = [0.,0.] #intersection point with bdry
    D = 0
    wall = None

    spec_case = False #bool for checking special case
    if x==xi and y==yi:
        spec_case = True
    if x==xi and y==yf:
        spec_case = True
    if x==xf and y==yi:
        spec_case = True
    if x==xf and y==yf:
        spec_case = True
    assert spec_case is False
    if x==xi: #starting on left boundry
        assert th<np.pi/2 or th>3*np.pi/2
        a1 = np.arctan((yf-y)/(xf-xi))
        a4 = np.arctan((xf-xi)/(y-yi)) + 3*np.pi/2
        if th<a1:
            P[0] = xf
            P[1] = np.tan(th)*(xf-x) + y
            wall = 'right'
        elif th<np.pi/2:
            P[1] = yf
            P[0] = (yf -y)/np.tan(th) + x
            wall = 'top'
        elif th<a4:
            P[1] = yi
            P[0] = (yi -y)/np.tan(th) + x
            wall = 'bottom'
        else:
            P[0] = xf
            P[1] = np.tan(th)*(xf-x) + y
            wall ='right'
        D = dist(p,P)
    elif x==xf: #right bdry
        assert th>np.pi/2 and th<3*np.pi/2
        a2 = np.arctan((xf-xi)/(yf-y)) + np.pi/2
        a3 = np.arctan(y/(xf-xi)) + np.pi
        if th<a2:
            P[1] = yf
            P[0] = (yf -y)/np.tan(th) + x
            wall ='top'
        elif th<a3:
            P[0] = xi
            P[1] = np.tan(th)*(xi-x) + y
            wall ='left'
        else:
            P[1] = yi
            P[0] = (yi -y)/np.tan(th) + x
            wall ='bottom'
        D = dist(p,P)
    elif y==yi: #top bdry
        assert th<np.pi
        a1 = np.arctan(yf/(xf-x))
        a2 = np.arctan(x/(yf-yi)) + np.pi/2
        # print(a1*180/np.pi)
        if th<a1:
            P[0] = xf
            P[1] = np.tan(th)*(xf-x) + y
            wall ='right'
        elif th<a2:
            P[1] = yf
            P[0] = (yf -y)/np.tan(th) + x
            wall ='top'
        else:
            P[0] = xi
            P[1] = np.tan(th)*(xi-x) + y
            wall = 'left'
        D = dist(p,P)
    elif y==yf: #bottom bdry
        assert th>np.pi
        a3 = np.arctan(yf/(x-xi)) + np.pi
        a4 = np.arctan((xf-x)/yf) + 3*np.pi/2
        if th<a3:
            P[0] = xi
            P[1] = np.tan(th)*(xi-x) + y
            wall = 'left'
        elif th<a4:
            P[1] = yi
            P[0] = (yi -y)/np.tan(th) + x
            wall ='bottom'
        else:
            P[0] = xf
            P[1] = np.tan(th)*(xf-x) + y
            wall ='right'
        D = dist(p,P)
    else: #inbetween
        a1 = np.arctan((yf-y)/(xf-x))
        a2 = np.arctan((x-xi)/(yf-y)) + np.pi/2
        a3 = np.arctan((y-yi)/(x-xi)) + np.pi
        a4 = np.arctan((xf-x)/(y-yi)) + 3*np.pi/2
        assert th!=a1 and th!=a2 and th!=a3 and th!=a4
        if th < a1 or th>a4:
            P[0] = xf
            P[1] = np.tan(th)*(xf-x) + y
            wall = 'right'
        elif th < a2:
            P[1] = yf
            P[0] = (yf-y)/np.tan(th) + x
            wall = 'top'
        elif th < a3:
            P[0] = xi
            P[1] = np.tan(th)*(xi-x) + y
            wall = 'left'
        else:
            P[1] = yi
            P[0] = (yi-y)/np.tan(th) + x
            wall = 'bottom'
        D = dist(p,P)
    D1 = 0 #consider if it deflects
    # print(wall)
    if not MT.from_bdry or len(MT.seg)>1: #if it's not from bdry or more than 1 segment
        # print('AA')
        D1 = MT.tip_l
    elif MT.from_bdry: #from bdry
        if rank == MT.origin: #if from same processor, can look at it previous MTs
            #TODO MIGHT NOT NEED THIS IF I REPLACE W/ prev_t
            idx = MT.prev_mt #index of prev MT
            l_idx = mt_to_l(MT_list, idx)#convert to list index
            MTp = MT_list[l_idx] #prev MT
        
            if MTp.from_bdry and len(MTp.seg) == 2: #if extension is extension
                tl = 0
                if MTp.origin == rank: #current MT is native to process
                    l_idx2 = mt_to_l(MT_list,MTp.prev_mt)#may need to consider ext of ext
                    MTp2 = MT_list[l_idx2]
                    tl = MTp.seg_dist[0] + (MTp2.update_t[-1] - MTp2.update_t[-2])#already grown
                    D1 = MT.tip_l - tl
                else: #previous MT is an extension from another process
                    tl = MT.update_t[-1] - MTp.prev_t #already grown
                    D1 = MT.tip_l - tl
                assert tl<= MT.tip_l
                assert D1>0
            else:
                D1 = MT.tip_l - (MTp.update_t[-1]- MTp.update_t[-2]) #tip length grown on new MT
                assert D1>0
        else:
            D1 = MT.tip_l - MT.update_t[-1]+MT.prev_t #length left to grow
            assert D1>0
    # D1 = 1
    if D1<D: #if deflects before intersection w/ bdry
        D = D1 #reassign
        wall = 'deflect' #keyword
        P[0] = p[0] + D1*np.cos(th) #deflected point
        P[1] = p[1] + D1*np.sin(th)
    if abs(P[0])>xdomain[1] or abs(P[1])>ydomain[1]:
        print('MT has gone crazy: ', MT.number)
        print('Point: ', P)
        assert abs(P[0])<=xdomain[1] and abs(P[1])<=ydomain[1]
    return(D,P,wall)

def compare(mt1,mt2,t):
    '''
    Parameters
    ----------
    mt1 : mt class. Check if m1 collides w/ m2 (ORDER MATTERS!)
    mt2 : mt class. To be checked if it's collided against
    t: Float. Current time
    Returns
    -------
    Float. shortest time of next collision
    '''
    output = compare_return() #initiate output
    tol = 1e-15 #numerical tolerance for intersection distance, otherwise no intersection
    if (mt1.hit_bdry is True and mt2.hit_bdry is True) or \
       (mt1.exist is False or mt2.exist is False) or (mt1.hit_bdry is True and mt2.grow is False) or\
       (mt2.hit_bdry is True and mt1.grow is False):
        output.dist = None
        output.policy = 'no_collision'
        output.point = None
        return(output)
    elif (mt1.grow is False and mt2.grow is False): #IF BOTH SHRINKING
        output.policy='no_collision'
    elif (mt1.grow and mt2.grow) is True: #FIRST ONLY LOOKING AT BOTH GROWING
        seg1, seg2 = mt1.seg, mt2.seg #assign segment points
        th1,th2 = mt1.angle, mt2.angle
        old_dist1, old_dist2 = mt1.seg_dist, mt2.seg_dist
        l1, l2 = len(seg1), len(seg2) #for indexing
        t_prev1, p1_prev = mt1.update_t[-1], seg1[l1-1] #last updated point and time
        t_prev2, p2_prev = mt2.update_t[-1], seg2[l2-1] #last updated point and time
        assert(len(old_dist1) == len(seg1)-1 and len(old_dist2) == len(seg2)-1)

        col_dist1t2 = [] #collision distances from 1 to 2
        col_dist2t1 = [] #coliision distances from 2 to 1

        point_1t2 = [] #store their respective collision locations
        point_2t1 = []

        #We now check if mt1 collides w/ any segments of mt2
        for i in range(l2):
            p2 = seg2[i] #point traj to be collided with
            col_result = inter(p1_prev,p2,th1[l1-1],th2[i]) #collision results
            if col_result[0] is True:
                d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                pt = col_result[1] #point of collision
                if i==l2-1: #checking collision w/ the two growing ends
                    d_g1 = d1 - (t-t_prev1) #distance at each mt collision
                    d_g2 = d2 - (t-t_prev2)
                    d_g = max(d_g1,d_g2) #distance grown by mt's is always max
                    if d_g >= tol: #cannot be an intersection that just happened
                        if d_g==d_g1: #store whichever collides with which
                            col_dist1t2.append(d_g)
                            point_1t2.append(pt)
                        else:
                            col_dist2t1.append(d_g)
                            point_2t1.append(pt)
                else:#check intersection of mt1 head w/ previous mt2 segments
                    if d2<= old_dist2[i]: #must be less than segment length for collision
                        d_g = d1-(t-t_prev1) #total distance grown since t
                        if d_g >= tol:
                            col_dist1t2.append(d_g) #store collision distance
                            point_1t2.append(pt)
        #check if m2 collides w/ any segments of mt1
        for j in range(l1):
            if j!=l1-1: #ignore collision w/ dynamic ends -already found above
                p1 = seg1[j] #point traj to be collided with
                col_result = inter(p2_prev,p1,th2[l2-1],th1[j]) #collision results
                pt = col_result[1]
                if col_result[0] is True:
                    d2, d1 = col_result[2][0], col_result[2][1] #distance to collision from mt2 and mt1
                    if d1<= old_dist1[j]: #must be less than segment length for collision
                        d_g = d2-(t-t_prev2) #total distance grown since t
                        if d_g >= tol:
                            col_dist2t1.append(d_g) #store collision distance
                            point_2t1.append(pt)
        if len(col_dist1t2) ==0 and len(col_dist2t1) ==0:
            output.policy = 'no_collision'
            return(output)
        elif len(col_dist1t2) ==0: #if no collisions
            which_min = 1
            i_2t1 = np.argmin(col_dist2t1)
        elif len(col_dist2t1) ==0:
            which_min = 0
            i_1t2 = np.argmin(col_dist1t2)
        else:
            min_1t2, i_1t2 = np.min(col_dist1t2), np.argmin(col_dist1t2) #find min distances
            min_2t1, i_2t1 = np.min(col_dist2t1), np.argmin(col_dist2t1)
            glob_min = [min_1t2,min_2t1] #array for total min
            which_min = np.argmin(glob_min) #find which occurs
        if which_min==0: #if mt1 hits mt2 first
            if col_dist1t2[i_1t2]<=0:
                output.policy = 'no_collision'
            else:
                output.policy = '1hit2'
                output.point = point_1t2[i_1t2]
                output.dist = col_dist1t2[i_1t2]
        else:
            if col_dist2t1[i_2t1]<=0:
                output.policy = 'no_collision'
            else:
                output.policy = '2hit1'
                output.point = point_2t1[i_2t1]
                output.dist = col_dist2t1[i_2t1]
    elif (not mt1.hit_bdry and not mt2.hit_bdry) is True: #IF ONE END IS SHRINKING
        assert mt1.grow is False or mt2.grow is False
        which = None #which tip hits
        if mt1.grow is False: #figure which is shrinking
            mts, mtg = mt1, mt2 #assign
            which = 2
        else:
            mts, mtg = mt2, mt1
            which = 1
        assert len(mts.seg) >= 2
        seg1, seg2 = mtg.seg, mts.seg #assign segment points
        th1,th2 = mtg.angle, mts.angle
        old_dist1, old_dist2 = mtg.seg_dist, mts.seg_dist
        l1, l2 = len(seg1), len(seg2) #for indexing
        t_prev1, p1_prev = mtg.update_t[-1], seg1[l1-1] #last updated point and time
        t_prev2, p2_prev = mts.update_t[-1], seg2[l2-1] #last updated point and time

        assert(len(old_dist1) == len(seg1)-1 and len(old_dist2) == len(seg2)-1)

        mt2_l = np.sum(old_dist2)#total length of shrinking mt, cannot got lower than this
        col_dist1t2 = [mt2_l/v_s - (t-t_prev2)] #collision distances from 1 to 2, growing can only collide w/ shrinking
        point_1t2 = [[0,0]]
        #We now check if mt1 collides w/ any segments of mt2
        for i in range(l2-1):
            if i==l2-2: #checking collision w/ the two dynamic ends
                p2 = seg2[i] #point traj to be collided with
                col_result = inter(p1_prev,p2,th1[l1-1],th2[i]) #collision results
                if col_result[0] is True:
                    d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                    if d2<= old_dist2[i]: #has to collide on segment
                        d_g = d1 - (t-t_prev1) #distance of mt1 collision, also time taken to grow this
                        d_s = v_s*(t-t_prev2+d_g) #distance shrank
                        d_segf = old_dist2[i] - d_s #total distance left on the segment
                        assert d2 != d_segf
                        if d2 < d_segf:
                            if d_g >= tol:
                                col_dist1t2.append(d_g)
                                point_1t2.append(col_result[1])
            else:#check intersection of mt1 head w/ previous mt2 segments
                p2 = seg2[i] #point traj to be collided with
                col_result = inter(p1_prev,p2,th1[l1-1],th2[i])
                if col_result[0] is True:
                    d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                    if d2<= old_dist2[i]: #has to collide on segment
                        d_g = d1 - (t-t_prev1) #distance/time taken to grow starting at t
                        d_s = v_s*(t-t_prev2+d_g) #distance shrank of mt2
                        d_segf = mt2_l - d_s #length of mt2 left
                        for k in range(i): #total length of mt2 at intersection point, add prev segs
                            d2 += old_dist2[k]
                        assert d2 != d_segf
                        if d2< d_segf: #if shrank to less than intersection distance
                            if d_g >= tol:
                                col_dist1t2.append(d_g)
                                point_1t2.append(col_result[1])
        if len(col_dist1t2)==1: #if not additional points are added, mt shrinks away before collision
            output.dist = col_dist1t2[0]
            if which == 1:
                output.policy = 'no_collision'#'2disap'
            else:
                output.policy = 'no_collision'#'1disap'
        else: #collision occurs
            min_1t2, i_1t2 = np.min(col_dist1t2), np.argmin(col_dist1t2)
            output.point = point_1t2[i_1t2]
            output.dist = col_dist1t2[i_1t2]
            if which == 2:
                output.policy = '2hit1'
            else:
                output.policy = '1hit2'
    elif (mt1.hit_bdry or mt2.hit_bdry) is True: #IF ONE IS GROWING AND ONE IS ON THE BORDER
        assert (mt1.grow or mt2.grow) is True
        which = None
        if mt1.hit_bdry is True: #figure which is growing
            mtb, mtg = mt1, mt2 #assign
            which = 2
        else:
            which = 1
            mtb, mtg = mt2, mt1
        assert len(mtb.seg) >= 2
        seg1, seg2 = mtg.seg, mtb.seg #assign segment points
        th1,th2 = mtg.angle, mtb.angle
        old_dist1, old_dist2 = mtg.seg_dist, mtb.seg_dist
        l1, l2 = len(seg1), len(seg2) #for indexing
        t_prev1, p1_prev = mtg.update_t[-1], seg1[l1-1] #last updated point and time
        t_prev2, p2_prev = mtb.update_t[-1], seg2[l2-1] #last updated point and time
        assert(len(old_dist1) == len(seg1)-1 and len(old_dist2) == len(seg2)-1)
        assert mt1.checkd() and mt2.checkd()
        col_dist1t2 = [] #collision distances from 1 to 2
        point_1t2 = [] #store their respective collision locations
        neg_dist = 0
        #We now check if mt1 collides w/ any segments of mt2
        for i in range(l2-1):
            p2 = seg2[i] #point traj to be collided with
            col_result = inter(p1_prev,p2,th1[l1-1],th2[i]) #collision results
            if col_result[0] is True:
                d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                pt = col_result[1]
                if d2<= old_dist2[i]: #must be less than segment length for collision
                    d_g = d1-(t-t_prev1) #total distance grown since t
                    if d_g<= tol:
                        neg_dist += 1
                    else:
                        col_dist1t2.append(d_g) #store collision distance
                        point_1t2.append(pt)
        if len(col_dist1t2)==0 or neg_dist>0:
            output.policy = 'no_collision'
        else:
            min_1t2, i_1t2 = np.min(col_dist1t2), np.argmin(col_dist1t2)
            output.point = point_1t2[i_1t2]
            output.dist = col_dist1t2[i_1t2]
            if which == 2:
                output.policy = '2hit1'
            else:
                output.policy = '1hit2'
    if output.policy != 'no_collision':
        assert output.dist > tol
            # output.policy = 'no_collision'
    return(output)

def which_seg(mt1,mt2,t):
    '''
    Given that mt1 collides with mt2, find which segment of mt2 it collides with
    Parameters
    ----------
    mt1 : mt class for MT that collides
    mt2 : mt class for MT that does not change

    Returns
    -------
    index of segment which get collided with

    '''
    index = None
    if (mt1.grow and mt2.grow) is True: #FIRST ONLY LOOKING AT BOTH GROWING
        seg1, seg2 = mt1.seg, mt2.seg #assign segment points
        th1,th2 = mt1.angle, mt2.angle
        old_dist1, old_dist2 = mt1.seg_dist, mt2.seg_dist
        l1, l2 = len(seg1), len(seg2) #for indexing
        t_prev1, p1_prev = mt1.update_t[-1], seg1[l1-1] #last updated point and time
        t_prev2, p2_prev = mt2.update_t[-1], seg2[l2-1] #last updated point and time
        assert(len(old_dist1) == len(seg1)-1 and len(old_dist2) == len(seg2)-1)

        col_dist1t2 = [] #collision distances from 1 to 2

        seg2_idx = []

        #We now check if mt1 collides w/ any segments of mt2
        for i in range(l2):
            p2 = seg2[i] #point traj to be collided with
            col_result = inter(p1_prev,p2,th1[l1-1],th2[i]) #collision results
            if col_result[0] is True:
                d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                if i==l2-1: #checking collision w/ the two growing ends
                    d_g1 = d1 - (t-t_prev1) #distance at each mt collision
                    d_g2 = d2 - (t-t_prev2)
                    d_g = max(d_g1,d_g2) #distance grown by mt's is always max
                    if d_g==d_g1: #we know mt1 intersects w/ mt2, don't care if mt2 intersects w/ mt1 here
                        col_dist1t2.append(d_g) #store distances and respective segment indices
                        seg2_idx.append(i)
                else:#check intersection of mt1 head w/ previous mt2 segments
                    if d2<= old_dist2[i]: #must be less than segment length for collision
                        d_g = d1-(t-t_prev1) #total distance grown since t
                        col_dist1t2.append(d_g) #store collision distance
                        seg2_idx.append(i)
        i_1t2 = np.argmin(col_dist1t2) #find min distances
        index = seg2_idx[i_1t2] #find which segment

    elif (not mt1.hit_bdry and not mt2.hit_bdry) is True: #IF ONE END IS SHRINKING
        assert mt1.grow is False or mt2.grow is False
        mts, mtg = mt2, mt1
        seg1, seg2 = mtg.seg, mts.seg #assign segment points
        th1,th2 = mtg.angle, mts.angle
        old_dist1, old_dist2 = mtg.seg_dist, mts.seg_dist
        l1, l2 = len(seg1), len(seg2) #for indexing
        t_prev1, p1_prev = mtg.update_t[-1], seg1[l1-1] #last updated point and time
        t_prev2 = mts.update_t[-1] #last updated point and time

        assert(len(old_dist1) == len(seg1)-1 and len(old_dist2) == len(seg2)-1)

        mt2_l = np.sum(old_dist2)#total length of shrinking mt, cannot got lower than this
        col_dist1t2 = [mt2_l/v_s - (t-t_prev2)] #collision distances from 1 to 2, growing can only collide w/ shrinking
        assert mt2_l/v_s - (t-t_prev2) >= 0
        seg2_idx = [0] #indices of collision segments
        #We now check if mt1 collides w/ any segments of mt2
        for i in range(l2-1):
            if i==l2-2: #checking collision w/ the two dynamic ends
                p2 = seg2[i] #point traj to be collided with
                col_result = inter(p1_prev,p2,th1[l1-1],th2[i]) #collision results
                if col_result[0] is True:
                    # print(col_result[1])
                    d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                    if d2<= old_dist2[i]: #has to collide on segment
                        d_g = d1 - (t-t_prev1) #distance of mt1 collision, also time taken to grow this
                        d_s = v_s*(t-t_prev2+d_g) #distance shrank
                        d_segf = old_dist2[i] - d_s #total distance left on the segment
                        assert d2 != d_segf
                        if d2 < d_segf:
                            col_dist1t2.append(d_g)
                            seg2_idx.append(i)
                            # point_1t2.append(col_result[1])
            else:#check intersection of mt1 head w/ previous mt2 segments
                p2 = seg2[i] #point traj to be collided with
                col_result = inter(p1_prev,p2,th1[l1-1],th2[i])
                if col_result[0] is True:
                    d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                    if d2<= old_dist2[i]: #has to collide on segment
                        d_g = d1 - (t-t_prev1) #distance/time taken to grow starting at t
                        d_s = v_s*(t-t_prev2+d_g) #distance shrank of mt2
                        d_segf = mt2_l - d_s #length of mt2 left
                        for k in range(i): #total length of mt2 at intersection point, add prev segs
                            d2 += old_dist2[k]
                        assert d2 != d_segf
                        if d2< d_segf: #if shrank to less than intersection distance
                            col_dist1t2.append(d_g)
                            seg2_idx.append(i)
        i_1t2 = np.argmin(col_dist1t2)
        index = seg2_idx[i_1t2]
    elif (mt1.hit_bdry or mt2.hit_bdry) is True: #IF ONE IS GROWING AND ONE IS ON THE BORDER
        assert (mt1.grow or mt2.grow) is True
        mtb, mtg = mt2, mt1
        assert len(mtb.seg) >= 2
        seg1, seg2 = mtg.seg, mtb.seg #assign segment points
        th1,th2 = mtg.angle, mtb.angle
        old_dist1, old_dist2 = mtg.seg_dist, mtb.seg_dist
        l1, l2 = len(seg1), len(seg2) #for indexing
        t_prev1, p1_prev = mtg.update_t[-1], seg1[l1-1] #last updated point and time
        t_prev2= mtb.update_t[-1] #last updated point and time

        assert(len(old_dist1) == len(seg1)-1 and len(old_dist2) == len(seg2)-1)
        # assert mt1.checkd() and mt2.checkd()
        col_dist1t2 = [] #collision distances from 1 to 2
        seg2_idx = []
        #We now check if mt1 collides w/ any segments of mt2
        for i in range(l2-1):
            p2 = seg2[i] #point traj to be collided with
            col_result = inter(p1_prev,p2,th1[l1-1],th2[i]) #collision results
            if col_result[0] is True:
                d1, d2 = col_result[2][0], col_result[2][1] #distance to collision from mt1 and mt2
                if d2<= old_dist2[i]: #must be less than segment length for collision
                    d_g = d1-(t-t_prev1) #total distance grown since t
                    col_dist1t2.append(d_g) #store collision distance
                    seg2_idx.append(i)
        i_1t2 = np.argmin(col_dist1t2)
        index = seg2_idx[i_1t2]
    return(index)

if __name__ == '__main__':
    p1 = [2,3]
    p2 = [0,0]
    th1, th2 = 3.5,1
    print(inter(p1,p2,th1,th2))
    M1 = mt(1)
    M1.seg = [[0.2,0.1]]
    # M1.seg_dist=[0.37989652999425766]
    M1.angle = [0.2]
    M1.tip_l = 200
    print(inter_bdry(M1,2))
