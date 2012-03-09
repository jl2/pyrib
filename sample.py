#!/usr/bin/env python3

from pyrib import *

import itertools
import math
import sys

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
        RiRotate(xrot, 1.0, 0.0, 0.0)
    else:
        RiRotate(-xrot, 1.0, 0.0, 0.0)

    if (direction[0] > 0):
        RiRotate(-yrot, 0.0, 1.0, 0.0)
    else:
        RiRotate(yrot, 0.0, 1.0, 0.0)


def PlaceCamera(cam):
    RiRotate(-cam['roll'], 0.0, 0.0, 1.0);
    la = cam['look_at']
    loc = cam['location']
    direction = [la[0]-loc[0],
                 la[1]-loc[1],
                 la[2]-loc[2]]
    AimZ(direction)
    negloc = [-1*x for x in loc]
    RiTranslate(*negloc)

def f(u,v):
    return 2*math.sin(u) * math.cos(v)

def main(args):
    rad = 20
    cam = {'location': [rad,rad,rad],
           'look_at': [0,0,0],
           'roll':0}
    
    RiBegin()
    RiSides(2);
    NUM_FRAMES = 20
    t = 0
    dt = 2*math.pi/(NUM_FRAMES-1);

    for i in range(NUM_FRAMES):
        RiFrameBegin(i);
        print("Rendering frame", i)
        RiDisplay("images/image{}.tif".format(i),"file","rgba")
        RiFormat(800, 600,  1.25)

        cam['location'] = [rad * math.sin(t), rad, rad * math.cos(t)]
        t += dt;

        RiProjection("perspective")
        PlaceCamera(cam)

        RiWorldBegin()
        umin = vmin = -4*math.pi
        umax = vmax = 4*math.pi
        steps = 128
        du = (umax - umin) / (steps-1)
        dv = (vmax - vmin) / (steps-1)

        RiColor(0,1,0)
        for pt in itertools.product(range(0,steps), range(0,steps)):
            u,v = pt
            u *= du
            v *= dv
            u += umin
            v += vmin
            RiTransformBegin()
            RiPolygon(3, P = [u,f(u,v),v,
                              u+du,f(u+du,v+dv),v+dv,
                              u+du,f(u+du,v),v,
                              ])
            RiPolygon(3, P = [u,f(u,v),v,
                              u,f(u,v+dv),v+dv,
                              u+du,f(u+du,v+dv),v+dv])
            RiTransformEnd()

        RiWorldEnd()
        RiFrameEnd()
    RiEnd();

if __name__=='__main__':
    main(sys.argv[1:])
