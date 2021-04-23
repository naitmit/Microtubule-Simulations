#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For determining the parameters of the simulation

@author: tim
"""
L = 80 #in micrometers - scaled s.t. this is 1
tip_length = 8 #in micrometers

tip_scaled = tip_length/L #scaled tip length
v_s = 0.5 #scaled speed, growth is set to 1

xdomain = [0,1] #scaled domain
ydomain = [0,1]
