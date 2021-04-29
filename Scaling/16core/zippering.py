#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 23:09:10 2021

@author: tim
"""
import numpy as np
import random as rnd

rnd.seed(1)
def zip_cat(angle1,angle2,pt,r):
    '''
    Determine whether the intersection results in catastrophe
    
    Parameters
    ----------
    angle1 : angle of tip which collides
    angle2 : angle of barrier MT
    pt : point of intersection

    Returns
    -------
    None.

    '''
    catastrophe = False
    th2,th1 = max(angle1,angle2),min(angle1,angle2)
    # print(angle1/np.pi,angle2/np.pi)
    d = 0# 2.5e-3 #distance away from MT for zippering
    th_crit = 2*np.pi/9 #critical angle
    new_pt = pt
    new_angle = angle1
    if th2 > 3*np.pi/2 and th1<np.pi/2:
        a1 = 2*np.pi-th2 #one angle
        a2 = th1
        b = a1+a2 #incident angle
        if th2 == angle1: #incoming angle is largest
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = np.pi + th1
                else: #catastrophe TODO
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
        else: 
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2-np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
    elif th2 > np.pi/2 and th2< np.pi and th1<np.pi/2:
        a1 = np.pi-th2 #one angle
        a2 = th1
        b = a1+a2 #incident angle
        if th2 == angle1: #incoming angle is largest
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1+np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
        else: 
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2+np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
    elif th2 > 3*np.pi/2 and th1<3*np.pi/2 and th1 > np.pi:
        a1 = th1- np.pi#one angle
        a2 = 2*np.pi - th2
        b = a1+a2 #incident angle
        if th2 == angle1: #incoming angle is largest
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1-np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
        else: 
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2-np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
    elif th2 > np.pi and th2<3*np.pi/2 and th1 < np.pi and th1 > np.pi/2:
        a1 = th2- np.pi#one angle
        a2 = np.pi-th1
        b = a1+a2 #incident angle
        if th2 == angle1: #incoming angle is largest
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1 + np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                    new_angle = th1
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
        else: 
            if b >= np.pi/2: #incident angle is large
                b2 = np.pi - b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2 - np.pi
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
            else: #incident angle is good
                b2 = b #redef incident angle
                if b2 <= th_crit: #zippering
                    new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                    new_angle = th2
                else: #catastrophe 
                    # r = rnd.randint(0, 1)
                    if r== 0:
                        catastrophe = True
    elif th2 < np.pi/2 and th1 < np.pi/2:
        b = th2-th1 #incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > np.pi and th2< 3*np.pi/2 and th1 < np.pi/2 and (th2-np.pi)>th1:
        a1 = np.pi-th2#one angle
        a2 = th1
        b = a1-a2#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1 + np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2 - np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > np.pi and th2< 3*np.pi/2 and th1 < np.pi/2 and (th2-np.pi)<th1:
        a1 = th1#one angle
        a2 = th2-np.pi
        b = a1-a2#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1 + np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2 - np.pi 
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > np.pi and th2< 3*np.pi/2 and th1 > np.pi and th1< 3*np.pi/2:
        b = th2-th1#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > np.pi/2 and th2< np.pi and th1 > np.pi/2 and th1 < np.pi:
        b = th2-th1#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1
            else: #catastrophe
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > 3*np.pi/2 and th1 > np.pi/2 and th1 < np.pi and (th1+np.pi)>th2:
        a1 = th1+np.pi#one angle
        a2 = th2
        b = a1-a2#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1 + np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2 - np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > 3*np.pi/2 and th1 > np.pi/2 and th1 < np.pi and (th1+np.pi)<th2:
        a1 = th2#one angle
        a2 = th1+np.pi
        b = a1-a2#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1 + np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2 - np.pi
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    elif th2 > 3*np.pi/2 and th1 > 3*np.pi/2:
        a1 = th2#one angle
        a2 = th1
        b = a1-a2#incident angle
        if th2 == angle1: #incoming angle is largest
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
                new_angle = th1
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
        else: 
            if b <= th_crit: #zippering
                new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
                new_angle = th2
            else: #catastrophe 
                # r = rnd.randint(0, 1)
                if r== 0:
                    catastrophe = True
    # else:
    #     print('ANGLES GONE WEIRD')
    #     print(th1/np.pi, th2/np.pi)
    # print('NEW ANGLE', new_angle/np.pi)
    return(new_angle, new_pt, catastrophe)

# def zip_cat(angle1,angle2,pt):
#     '''
#     Parameters
#     ----------
#     angle1 : angle of tip which collides
#     angle2 : angle of barrier MT
#     pt : point of intersection

#     Returns
#     -------
#     None.

#     '''
#     th2,th1 = max(angle1,angle2),min(angle1,angle2)
#     # print(angle1/np.pi,angle2/np.pi)
#     d = 0 #distance away from MT for zippering
#     th_crit = 2*np.pi/9 #critical angle
#     new_pt = pt
#     new_angle = angle1
#     if th2 > 3*np.pi/2 and th1<np.pi/2:
#         a1 = 2*np.pi-th2 #one angle
#         a2 = th1
#         b = a1+a2 #incident angle
#         if th2 == angle1: #incoming angle is largest
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = np.pi + th1
#                 # else: #catastrophe TODO
#                 #     r = rnd.randint(0, 1)
#                 #     if r== 0:
#                 #         new_angle = 
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1
#         else: 
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2-np.pi
#                 # else: #catastrophe TODO
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2
#     elif th2 > np.pi/2 and th2< np.pi and th1<np.pi/2:
#         a1 = np.pi-th2 #one angle
#         a2 = th1
#         b = a1+a2 #incident angle
#         if th2 == angle1: #incoming angle is largest
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1
#                 # else: #catastrophe TODO
#                 #     r = rnd.randint(0, 1)
#                 #     if r== 0:
#                 #         new_angle = 
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1+np.pi
#         else: 
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2
#                 # else: #catastrophe TODO
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2+np.pi
#     elif th2 > 3*np.pi/2 and th1<3*np.pi/2 and th1 > np.pi:
#         a1 = th1- np.pi#one angle
#         a2 = 2*np.pi - th2
#         b = a1+a2 #incident angle
#         if th2 == angle1: #incoming angle is largest
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1
#                 # else: #catastrophe TODO
#                 #     r = rnd.randint(0, 1)
#                 #     if r== 0:
#                 #         new_angle = 
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1-np.pi
#         else: 
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2
#                 # else: #catastrophe TODO
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2-np.pi
#     elif th2 > np.pi and th2<3*np.pi/2 and th1 < np.pi and th1 > np.pi/2:
#         a1 = th2- np.pi#one angle
#         a2 = np.pi-th1
#         b = a1+a2 #incident angle
#         if th2 == angle1: #incoming angle is largest
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1 + np.pi
#                 # else: #catastrophe TODO
#                 #     r = rnd.randint(0, 1)
#                 #     if r== 0:
#                 #         new_angle = 
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th2),pt[1] - d*np.sin(th2)] #new point
#                     new_angle = th1
#         else: 
#             if b >= np.pi/2: #incident angle is large
#                 b2 = np.pi - b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2 - np.pi
#                 # else: #catastrophe TODO
#             else: #incident angle is good
#                 b2 = b #redef incident angle
#                 if b2 <= th_crit: #zippering
#                     new_pt = [pt[0] - d*np.cos(th1),pt[1] - d*np.sin(th1)] #new point
#                     new_angle = th2
#     # else:
#     #     print('ANGLES GONE WEIRD')
#     # print('NEW ANGLE', new_angle/np.pi)
#     return(new_angle, new_pt)