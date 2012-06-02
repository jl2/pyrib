#!/usr/bin/env python3

# sample.py

# Copyright (c) 2012, Jeremiah LaRocco jeremiah.larocco@gmail.com

# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.

# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


import sys
sys.path.append('..')

import ri

import itertools
import math

def AimZ(direction):

    if (direction[0]==0 and direction[1]==0 and direction[2]==0):
        return

    xzlen = math.sqrt(direction[0]*direction[0]+direction[2]*direction[2])
    if (xzlen == 0):
        if (direction[1] < 0):
            yrot = 180
        else:
            yrot = 0
    else:
        yrot = 180.0*math.acos(direction[2]/xzlen)/math.pi

    yzlen = math.sqrt(direction[1]*direction[1]+xzlen*xzlen)
    xrot = 180*math.acos(xzlen/yzlen)/math.pi

    if (direction[1] > 0):
        ri.Rotate(xrot, 1.0, 0.0, 0.0)
    else:
        ri.Rotate(-xrot, 1.0, 0.0, 0.0)

    if (direction[0] > 0):
        ri.Rotate(-yrot, 0.0, 1.0, 0.0)
    else:
        ri.Rotate(yrot, 0.0, 1.0, 0.0)


def PlaceCamera(cam):
    ri.Rotate(-cam['roll'], 0.0, 0.0, 1.0);
    la = cam['look_at']
    loc = cam['location']
    direction = [la[0]-loc[0],
                 la[1]-loc[1],
                 la[2]-loc[2]]
    AimZ(direction)
    negloc = [-1*x for x in loc]
    ri.Translate(*negloc)

def f(u,v):
    return 2*math.sin(u) * math.cos(v)

def main(args):
    rad = 20
    cam = {'location': [rad,rad,rad],
           'look_at': [0,0,0],
           'roll':0}
    
    ri.Begin()
    ri.Sides(2);
    NUM_FRAMES = 20
    t = 0
    dt = 2*math.pi/(NUM_FRAMES-1);

    for i in range(NUM_FRAMES):
        ri.FrameBegin(i);
        print("Rendering frame", i)
        ri.Display("images/image{}.tif".format(i),"file","rgba")
        ri.Format(800, 600,  1.25)

        cam['location'] = [rad * math.sin(t), rad, rad * math.cos(t)]
        t += dt;

        ri.Projection("perspective")
        PlaceCamera(cam)

        ri.WorldBegin()
        umin = vmin = -4*math.pi
        umax = vmax = 4*math.pi
        steps = 128
        du = (umax - umin) / (steps-1)
        dv = (vmax - vmin) / (steps-1)

        ri.Color(0,1,0)
        for pt in itertools.product(range(0,steps), range(0,steps)):
            u,v = pt
            u *= du
            v *= dv
            u += umin
            v += vmin
            ri.TransformBegin()
            ri.Polygon(3, P = [u,f(u,v),v,
                              u+du,f(u+du,v+dv),v+dv,
                              u+du,f(u+du,v),v,
                              ])
            ri.Polygon(3, P = [u,f(u,v),v,
                              u,f(u,v+dv),v+dv,
                              u+du,f(u+du,v+dv),v+dv])
            ri.TransformEnd()

        ri.WorldEnd()
        ri.FrameEnd()
    ri.End();

if __name__=='__main__':
    main(sys.argv[1:])
