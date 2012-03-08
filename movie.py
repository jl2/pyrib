#!/usr/bin/env python3

from pyrib import *


import itertools
import argparse
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
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--frames', help='The number of frames',
                        type=int, default=24)
    parser.add_argument('--prefix', help="The base name of the images.  Will look for: \"images/basename*.tif\"",
                        default="sample")
    parser.add_argument('--out', help='The name of the output file.  Must be recognized by ffmpeg',
                        type=os.path.abspath, required=True, metavar='<outfile>')
    parser.add_argument('--audio', help='The audio file for the sound.  Must be supported by ffmpeg',
                        type=os.path.abspath, default=None, metavar='<audio file>')
    pargs = parser.parse_args(args)    

    rad = 16
    cam = {'location': [rad,rad,rad],
           'look_at': [0,0,0],
           'roll':0}
    
    RiBegin()
    RiSides(2);
    NUM_FRAMES = pargs.frames
    t = 0
    dt = 2*math.pi/(NUM_FRAMES-1);

    for i in range(NUM_FRAMES):
        RiFrameBegin(i);
        print("Rendering frame", i)
        RiDisplay("images/{}{:05d}.jpg".format(pargs.prefix, i),"jpeg","rgb")
        RiFormat(800, 600,  1.25)

        cam['location'] = [rad * math.sin(t), rad, rad * math.cos(t)]
        t += dt;

        RiProjection("perspective")
        PlaceCamera(cam)
        RiAttribute("light", "string shadow", "on")
        RiLightSource("distantlight", "from", to_fa([-40,80,40]), "to", to_fa([0,0,0]))
        RiImager("background")
        RiWorldBegin()
        umin = vmin = -4*math.pi
        umax = vmax = 4*math.pi
        steps = 128
        du = (umax - umin) / (steps-1)
        dv = (vmax - vmin) / (steps-1)

        RiAttributeBegin()
        RiColor(0,1,1)
        RiSurface("plastic" )
        RiTorus(8, 1, 0, 360, 360)
        RiAttributeEnd()
        
        RiAttributeBegin()
        RiColor(0,1,0)
        RiSurface("matte")
        RiParaboloid(4, 0, 4, 360)
        RiAttributeEnd()
        
        RiAttributeBegin()
        RiColor(0,1,0)
        RiSurface("plastic")
        RiHyperboloid([4,4,4], [4,4, 0], 360)
        RiAttributeEnd()
        
        RiAttributeBegin()
        RiColor(1,0,0)
        RiTranslate(0, 0, 8)
        RiSurface("checker")
        RiSphere(3, -3, 3, 360)
        RiAttributeEnd()
        
        RiWorldEnd()
        RiFrameEnd()
    RiEnd();

    out_movie = pargs.out
    idx = out_movie.rfind('.')
    tmp_movie = out_movie[:idx] + '_tmp' + out_movie[idx:]
    movie_cmd = 'ffmpeg -r 30 -b 2400 -i {directory}/{prefix}%05d.jpg -q 2 {outfile}'.format(directory='images',
                                                                                                 prefix=pargs.prefix,
                                                                                                 outfile=tmp_movie)
    print('Executing: {}'.format(movie_cmd))
    os.system(movie_cmd)
    if 'flac' in out_movie[idx:]:
        out_codec = 'mp2'
    else:
        out_codec = 'copy'

    if pargs.audio is not None:
        snd_cmd = 'ffmpeg -i {inmovie} -i "{insound}" -acodec {codec} -vcodec copy -ab 320k -ar 44100 {outfile}'.format(inmovie=tmp_movie,
                                                                                                                        insound=pargs.audio,
                                                                                                                        outfile=out_movie,
                                                                                                                        codec=out_codec)
        print('Executing: {}'.format(snd_cmd))
        os.system(snd_cmd)
    else:
        os.rename(tmp_movie, out_movie)

if __name__=='__main__':
    main(sys.argv[1:])
