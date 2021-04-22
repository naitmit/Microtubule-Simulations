#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 11:41:58 2021

@author: tim
"""
import matplotlib.pyplot as plt
import numpy as np
from parameters import v_s, xdomain, ydomain

def plot_snap(mt_list,t,k,dest='./plots3/',save=True):
    l = len(mt_list)
    for i in range(l):
        coord = mt_list[i].seg #get coordinates of segments
        x_seg = [p[0] for p in coord]
        y_seg = [p[1] for p in coord]
        th = mt_list[i].angle[-1] #angle
        t_i = mt_list[i].update_t[-1] #last update time
        if mt_list[i].exist: #if it exists, plot segments
            plt.scatter(x_seg,y_seg,s=2, color='black') #scatter of all segment points
        if t_i == t: #if it's been updated at this time
            if mt_list[i].exist:
                if mt_list[i].grow or mt_list[i].hit_bdry: #growing or stationary
                    x_i = x_seg[:-1] #all but most recent update
                    y_i = y_seg[:-1]
                    x_f = x_seg[-2:] #updated segment
                    y_f = y_seg[-2:]
                    plt.plot(x_i,y_i,color='green',linewidth=0.5)
                    if len(x_f)>1:
                        plt.plot(x_f,y_f,color='red',linewidth=3, label='MT event')
                else:
                    plt.plot(x_seg,y_seg,color='green',linewidth=0.5) #started to shrink
            else:
                plt.scatter([x_seg[0]],[y_seg[0]], s=10,color='purple') #yellow points are dissapeared MTs
        else: #if it's continuing to grow at this point in time
            if mt_list[i].exist:
                if mt_list[i].grow or mt_list[i].hit_bdry:
                    x2 = x_seg[-1] + np.cos(th)*(t-t_i)
                    y2 = y_seg[-1] + np.sin(th)*(t-t_i)
                    if x2 > xdomain[1]: #for plotting the tip, so that it's within the domain
                        x2 = xdomain[1]
                    elif x2<xdomain[0]:
                        x2 = xdomain[0]
                    
                    if y2 > ydomain[1]:
                        y2 = ydomain[1]
                    elif y2<ydomain[0]:
                        y2 = ydomain[0]
                        
                    x_f = [x_seg[-1], x2] #where it has grown to
                    y_f = [y_seg[-1], y2]
                    plt.plot(x_seg,y_seg,color='green',linewidth=0.5)
                    plt.plot(x_f,y_f, color='blue', linewidth=0.5)
                else:
                    assert not mt_list[i].grow
                    cumsum = np.append([0],np.cumsum(mt_list[i].seg_dist)) #cumulative sums
                    d = (t- t_i)*v_s #distance traversed
                    left_over = cumsum[-1] - d
                    j = 0
                    while left_over > cumsum[j]: #find index within cumsum
                        j+=1
                        if j > len(cumsum)-1:
                            break
                    j -= 1
                    delta = left_over - cumsum[j] #component sticking out from segment
                    x_i = x_seg[:j+1] #untouched points in shrinkage
                    y_i = y_seg[:j+1]
                    th = mt_list[i].angle[j]
                    x_f = [x_seg[j], x_seg[j] + np.cos(th)*delta] #where it has grown to
                    y_f = [y_seg[j], y_seg[j] + np.sin(th)*delta]
                    plt.plot(x_i,y_i,color='green',linewidth=0.5)
                    plt.plot(x_f,y_f,color='cyan',linewidth=.5)
            else:
                plt.scatter([x_seg[0]],[y_seg[0]],s=10,color='purple')
                # print('sdsd')
                
    plt.xlim(xdomain)
    plt.ylim(ydomain)
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