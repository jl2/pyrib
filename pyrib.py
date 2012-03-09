#!/usr/bin/env python3

# pyrib.py

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

# This is a thin wrapper around the Renderman library using ctypes
# It needs a lot of work

import ctypes
import sys
import os

delight_path = os.getenv("DELIGHT")
libnames = {'win32': '/bin/3delight.dll',
            'cygwin': '/bin/3delight.dll',
            'linux2': '/lib/lib3delight.so',
            'darwin': '/lib/lib3delight.dylib',
            # Try something...
            'default': '/lib/lib3delight.so',}

if sys.platform in libnames.keys():
    delight_lib = delight_path + libnames[sys.platform]
else:
    delight_lib = delight_path + libnames['default']

rlib = ctypes.CDLL(delight_lib)

RI_FRAMEBUFFER = 'framebuffer'
RI_FILE = 'file'
RI_ZFILE = 'zfile'

RI_RGB = 'rgb'
RI_RGBA = 'rgba'
RI_RGBZ = 'rgbz'
RI_RGBAZ = 'rgbaz'

RI_A = 'a'
RI_Z = 'z'
RI_AZ = 'az'

RI_PERSPECTIVE = 'perspective'
RI_ORTHOGRAPHIC = 'orthographic'

RI_HIDDEN = 'hidden'
RI_PAINT = 'paint'
RI_RAYTRACE = 'raytrace'
RI_PHOTON = 'photon'

RI_CONSTANT = 'constant'
RI_SMOOTH = 'smooth'

RI_FLATNESS = 'flatness'
RI_FOV = 'fov'

RI_AMBIENTLIGHT = 'ambientlight'
RI_POINTLIGHT = 'pointlight'
RI_DISTANTLIGHT = 'distantlight'
RI_SPOTLIGHT = 'spotlight'

RI_INTENSITY = 'intensity'
RI_LIGHTCOLOR = 'lightcolor'
RI_FROM = 'from'
RI_TO = 'to'
RI_CONEANGLE = 'coneangle'
RI_CONEDELTAANGLE = 'conedeltaangle'
RI_BEAMDISTRIBUTION = 'beamdistribution'

RI_MATTE = 'matte'
RI_METAL = 'metal'
RI_SHINYMETAL = 'shinymetal'
RI_PLASTIC = 'plastic'
RI_PAINTEDPLASTIC = 'paintedplastic'

RI_KA = 'Ka'
RI_KD = 'Kd'
RI_KR = 'Kr'
RI_KS = 'Ks'
RI_ROUGHNESS = 'roughness'

RI_TEXTURENAME = 'texturename'
RI_SPECULARCOLOR = 'specularcolor'

RI_DEPTHCUE = 'depthcue'
RI_FOG = 'fog'
RI_BUMPY = 'bumpy'

RI_MINDISTANCE = 'mindistance'
RI_MAXDISTANCE = 'maxdistance'
RI_BACKGROUND = 'background'
RI_DISTANCE = 'distance'
RI_AMPLITUDE = 'amplitude'

RI_RASTER = 'raster'
RI_SCREEN = 'screen'
RI_CAMERA = 'camera'
RI_WORLD = 'world'
RI_OBJECT = 'object'

RI_INSIDE = 'inside'
RI_OUTSIDE = 'outside'
RI_LH = 'lh'
RI_RH = 'rh'

RI_P = 'P'
RI_PZ = 'Pz'
RI_PW = 'Pw'
RI_N = 'N'
RI_NP = 'Np'
RI_CS = 'Cs'
RI_OS = 'Os'
RI_S = 's'
RI_T = 't'
RI_ST = 'st'

RI_BILINEAR = 'bilinear'
RI_BICUBIC = 'bicubic'

RI_PRIMITIVE = 'primitive'
RI_INTERSECTION = 'intersection'
RI_UNION = 'union'
RI_DIFFERENCE = 'difference'

RI_PERIODIC = 'periodic'
RI_NONPERIODIC = 'nonperiodic'
RI_CLAMP = 'clamp'
RI_BLACK = 'black'

RI_IGNORE = 'ignore'
RI_PRINT = 'print'
RI_ABORT = 'abort'
RI_HANDLER = 'handler'

RI_BOUNDS = 'bounds'

RI_LIMITS = 'limits'
RI_SHADOW = 'shadow'
RI_BIAS0 = 'bias0'
RI_BIAS1 = 'bias1'
RI_SAMPLE = 'sample'

RI_SEARCHPATH = 'searchpath'
RI_SHADER = 'shader'
RI_TEXTURE = 'texture'
RI_DISPLAY = 'display'

RI_WIDTH = 'width'
RI_CONSTANTWIDTH = 'constantwidth'

RI_COMMENT = 'comment'
RI_STRUCTURE = 'structure'
RI_VERBATIM = 'verbatim'

RI_IDENTIFIER = 'identifier'
RI_NAME = 'name'
RI_SHADINGGROUP = 'shadinggroup'

RI_SEARCHPATH = 'searchpath'
RI_SHADER = 'shader'
RI_TEXTURE = 'texture'


RI_LINEAR = 'linear'
RI_CUBIC = 'cubic'


RI_BUCKETSENDER = 'bucketsender'
RI_FUNCTION = 'function'

RI_THRESHOLD = 'threshold'
RI___THRESHOLD = '__threshold'

RI_HANDLEID = 'handleid'
RI___HANDLEID = '__HANDLEID'


RI_QUANTIZE = 'QUANTIZE'
RI_DITHER = 'dither'
RI_EXPOSURE = 'EXPOSURE'
RI_FILTER = 'filter'
RI_FILTERWIDTH = 'filterwidth'


RiBezierBasis = rlib.RiBezierBasis
RiBSplineBasis = rlib.RiBezierBasis
RiCatmullRomBasis = rlib.RiBezierBasis
RiHermiteBasis = rlib.RiBezierBasis
RiPowerBasis = rlib.RiBezierBasis

RI_BEZIERSTEP = 3
RI_BSPLINESTEP = 1
RI_CATMULLROMSTEP = 1
RI_HERMITESTEP = 2
RI_POWERSTEP = 4


rf = ctypes.c_float

RtColor = rf*3
RtPoint = rf*3
RtVector = rf*3
RtNormal = rf*3
RtHpoint = rf*4
RtBound = rf*6
RtInterval = rf*2

# typedef RtFloat RtMatrix[4][4];
# typedef RtFloat RtBasis[4][4];
def to_fa(val):
    return list(map(float, val))

def to_ia(val):
    return list(map(int, val))

def to_rs(val):
    return ctypes.c_char_p(val.encode('utf8'))

def rpa(val):
    rv = (RtPoint*len(val))(*list(map(lambda x: RtPoint(*x), val)))
    return rv

def to_rfa(val):
    rv = (ctypes.c_float*len(val))(*val)
    return rv

def to_ria(val):
    return (ctypes.c_int*len(val))(*list(map(ctypes.c_int, val)))

def to_rsa(val):
    return list(map(to_rs, val))

def to_rf(val):
    return ctypes.c_float(val)

def to_ri(val):
    return ctypes.c_int(val)

# This is super ugly, but works
# Most renderman functions take a variable number of arguments using "tokens" and "values", followed
# by RI_NULL.  It looks like this:
# RtFloat verts[] = {-10,-10,0,
#                     10,-10,0,
#                     10,10,0,
#                     -10,10,0};
# RiPolygon(4, "P", verts, RI_NULL)

# For each function there's a corresponding *V function that takes the an integer and the tokens and parameters separated into
# arrays.  The corresponding function looks like this:
# RtFloat verts[] = {-10,-10,0,
#                     10,-10,0,
#                     10,10,0,
#                     -10,10,0};
# RtToken tokens[] = {"P"};
# RtPointer parms[] = {verts};
# RiPolygonV(4, 1, tokens, parms);

# The following functions (to_c, get_c_type, make_void_p and parse_extra_args) convert a Python **kwargs parameter
# into the RtInt n, RtToken tokens[], RtPointer parms[] parameters used by the *V functions

# to_c() may be useful by itself

def to_c(val):
    tmp_type = type(val)
    if tmp_type == type(1):
        return ctypes.c_int(val)
    if tmp_type == type(1.2):
        return ctypes.c_float(val)
    if tmp_type == type('abc'):
        return ctypes.c_char_p(val.encode('utf8'))
    if tmp_type == type(list()):
        if len(val) >0:
            return (get_c_type(val[0])*len(val))(*list(map(to_c, val)))
        return None
    return val

def get_c_type(val):
    tmp_type = type(val)
    if tmp_type == type(1):
        return ctypes.c_int
    if tmp_type == type(1.2):
        return ctypes.c_float
    if tmp_type == type('abc'):
        return ctypes.c_char_p
    if tmp_type == type(list()):
        if len(val) >0:
            return get_c_type(val[0])*len(val)
        return None

def make_void_p(val):
    return ctypes.cast(to_c(val), ctypes.c_void_p)

def to_rfa_flatten(val):
    rv = []
    for v in val:
        if type(v) == type(list()):
            rv = rv + to_rfa(v)
        else:
            rv.append(float(v))
    return list(map(ctypes.c_float, rv))

# This function attempts to convert parameters based on expected types
def convert(token, val):
    if token not in {"P", "Pz", "Pw", "N", "Np", "Cs", "Os", "s", "t", "st"}:
        return to_c(val)
    
    rv = to_rfa_flatten(val)
    if token in {"P", "N", "Np", "Os", "Cs"} and len(rv)%3 != 0:
            raise RuntimeError('"{}" requires a multiple of three arguments!')

    if token in {"st"} and len(rv)%2 != 0:
        raise RuntimeError('"{}" requires a multiple of two arguments!')
    return (ctypes.c_float*len(rv))(*rv)

def parse_extra_args(args, kwargs):
    toks = []
    parms = []
    num_params = len(args) + len(kwargs.keys())
    if num_params == 0:
        return (0,None, None)

    for i in range(len(args)//2):
        toks.append(tok2c(args[i*2]))
        parms.append(ctypes.cast(to_c(args[i*2+1]), ctypes.c_void_p))

    
    i = 0
    for k in kwargs.keys():
        toks.append(tok2c(k))
        rib_val = convert(k, kwargs[k])
        parms.append(ctypes.cast(rib_val, ctypes.c_void_p))
        i += 1
    toks = (ctypes.c_char_p*len(toks))(*toks)
    parms = (ctypes.c_void_p*len(parms))(*parms)
    rv = (to_ri(len(toks)), toks, parms)
    return rv

def tok2c(val):
    return ctypes.c_char_p(val.encode('utf8'))

def RiBegin(name=''):
    if len(name)==0:
        rlib.RiBegin(0)
    else:
        rlib.RiBegin(tok2c(cname))

def RiEnd():
    rlib.RiEnd()


def RiCamera(camera):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiCameraV(tok2c(camera), n, toks, parms)

# DL_INTERFACE RtVoid RiCameraV(
# 		RtToken camera,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

def RiGetContext():
    return rlib.RiGetContext()
    

def RiContext(handle):
    rlib.RiContext(handle)

def RiFrameBegin(frame):
    rlib.RiFrameBegin(ctypes.c_int(frame))

def RiFrameEnd():
    rlib.RiFrameEnd()


def RiMotionBegin(n, *args):
    RiMotionBeginV(n, args)

def RiMotionBeginV(n, times):
    c_times_type = RtFloat * len(times)
    rlib.RiMotionBeginV(n, c_times_type(*times))

def RiMotionEnd():
    rlib.RiMotionEnd()

def RiSolidBegin(operation):
    rlib.RiSolidBegin(tok2c(operation))

def RiSolidEnd():
    rlib.RiSolidEnd()

def RiWorldBegin():
    rlib.RiWorldBegin()

def RiWorldEnd():
    rlib.RiWorldEnd()

def RiObjectBegin(*args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    return rlib.RiObjectBeginV(n, toks, parms)

def RiObjectEnd():
    rlib.RiObjectEnd()


def RiObjectInstance(handle):
    rlib.RiObjectInstance(handle)

def RiResource(handle, rtype, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)    
    rlib.RiResourceV(tok2c(handle), tok2c(rtype), n, toks, parms)

def RiResourceBegin():
    rlib.RiResourceBegin()

def RiResourceEnd():
    rlib.RiResourceEnd()

def RiFormat(xres, yres, aspect):
    rlib.RiFormat(ctypes.c_int(xres), ctypes.c_int(yres), to_rf(aspect))
    
def RiFrameAspectRatio(aspect):
    rlib.RiFrameAspectRatio(to_rf(aspect))

def RiScreenWindow(left, right, bottom, top):
    rlib.RiScreenWindow(to_rf(left), to_rf(right), to_rf(bottom), to_rf(top))

def RiClipping(hither, yon):
    rlib.RiClipping(to_rf(hither), to_rf(yon))

def RiClippingPlane(x, y, z, nx, ny, nz ):
    rlib.RiClippingPlane(to_rf(x), to_rf(y), to_rf(z), to_rf(nx), to_rf(ny), to_rf(nz))

def RiCropWindow(xmin, xmax, ymin, ymax):
    rlib.RiCropWindow(to_rf(xmin), to_rf(xmax), to_rf(ymin), to_rf(ymax))
    
def RiDepthOfField(fstop, focallength, focaldistance):
    rlib.RiDepthOfField(to_rf(fstop), to_rf(focallength), to_rf(focaldistance))

def RiProjection(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)    
    rlib.RiProjectionV(tok2c(name), n, toks,parms)

def RiShutter(smin, smax):
    rlib.RiShutter(to_rf(smin), to_rf(smax))

def RiDisplay(name, dtype, mode, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiDisplayV(tok2c(name), tok2c(dtype), tok2c(mode), n, toks, parms)
                  
def RiDisplayChannel(channel, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiDisplayChannelV(tok2c(channel), n,toks, parms)

def RiExposure(gain, gamma):
    rlib.RiExposure(to_rf(gain), to_rf(gamma))

def RiImager(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiImagerV(tok2c(name), n,toks, parms)

# TODO: Implement RiPixelFilter
# DL_INTERFACE RtVoid RiPixelFilter(RtFilterFunc filterfunc,
# 							   RtFloat xwidth, RtFloat ywidth);
def RiPixelSamples(xsamples, ysamples):
    rlib.RiPixelSamples(to_rf(xsamples), to_rf(ysamples))

def RiPixelVariance(variation):
    rlib.RiPixelVariance(to_rf(variation))

def RiQuantize(qtype, one, qmin, qmax, ampl):
    rlib.RiQuantize(tok2c(qtype), ctypes.c_int(one), ctypes.c_int(qmin), ctypes.c_int(qmax), to_rf(ampl));

def RiConcatTransform(transform):
     rlib.RiConcatTransform(to_c(transform))

def RiCoordinateSystem(space):
    rlib.RiCoordinateSystem(tok2c(space))

def RiScopedCoordinateSystem(space):
    rlib.RiScopedCoordinateSystem(tok2c(space))

def RiCoordSysTransform(space):
    rlib.RiCoordSysTransform(tok2c(space))

def RiIdentity():
    rlib.RiIdentity()

def RiPerspective(fov):
    rlib.RiPerspective(ctypes.c_float(fov))

def RiRotate(angle, dx, dy, dz):
    rlib.RiRotate(to_rf(angle), to_rf(dx), to_rf(dy), to_rf(dz))
def RiScale(dx, dy, dz):
    rlib.RiScale(to_rf(dx), to_rf(dy), to_rf(dz))
def RiTranslate(dx, dy, dz):
    rlib.RiTranslate(to_rf(dx), to_rf(dy), to_rf(dz))

def RiSkew(angle, dx1, dy1, dz1, dx2, dy2, dz2):
    rlib.RiSkew(to_rf(angle), to_rf(dx1), to_rf(dy1), to_rf(dz1),
                to_rf(dx2), to_rf(dy2), to_rf(dz2))

def RiTransform(transform):
    rlib.RiTransform(to_c(transform))

def RiTransformBegin():
    rlib.RiTransformBegin()

def RiTransformEnd():
    rlib.RiTransformEnd()

# TODO: Make sure this works
# DL_INTERFACE RtPoint* RiTransformPoints(RtToken fromspace,
# 									 RtToken tospace,
# 									 RtInt n, RtPoint points[]);
def RiTransformPoints(fromspace, tospace, n, points):
    return rlib.RiTransformPoints(tok2c(fromspace), tok2c(tospace), to_ri(n), to_c(points))

def RiAtmosphere(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiAtmosphereV(tok2c(name), n, toks, parms)

def RiDeformation(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiDeformationV(tok2c(name), n, toks, parms)

def RiDisplacement(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiDisplacementV(tok2c(name), n,toks,parms)

def RiExterior(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiExteriorV(tok2c(name), n, toks, parms)

def RiIlluminate(light, onoff):
    rlib.RiIlluminate(light, ctypes.c_bool(onoff))

def RiInterior(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiInteriorV(tok2c(name), n, toks, parms)

def RiShader(name, handle, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiShaderV(tok2c(name), to_c(handle), n, toks, parms)

def RiMatte(onoff):
    rlib.RiMatte(ctypes.c_bool(onoff))

def RiMultiplyShadingRate(ratemultiplier):
    rlib.RiMultiplyShadingRate(to_rf(ratemultiplier))

def RiShadingRate(size):
    rlib.RiShadingRate(to_rf(size))

def RiShadingInterpolation(itype):
    rlib.RiShadingInterpolation(tok2c(itype))

def RiSurface(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiSurfaceV(tok2c(name), n, toks, parms)

# TODO: implement RiArchiveRecord
# def RiArchiveRecord(RtToken type, __RI_CONST char *format, ...):
#     rlib.RiArchiveRecord(RtToken type, __RI_CONST char *format, ...)

# TODO: implement RiReadArchive
# DL_INTERFACE RtVoid RiReadArchive(
# 		RtToken name, RtArchiveCallback callback, ...);
# DL_INTERFACE RtVoid RiReadArchiveV(
# 		RtToken name, RtArchiveCallback callback,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

def RiArchiveBegin(archivename, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    return rlib.RiArchiveBegin(tok2c(archivename), n, toks, parms)

def RiArchiveEnd():
    rlib.RiArchiveEnd()


def RiIfBegin(expression, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiIfBeginV(tok2c(expression), n, toks, parms)

def RiElseIf(expression, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiElseIfV(tok2c(expression), n, toks, parms)

def RiElse():
    rlib.RiElse()

def RiIfEnd():
    rlib.RiIfEnd()

def RiAttribute(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiAttributeV(to_c(name), n, toks, parms)

def RiAttributeBegin():
    rlib.RiAttributeBegin()

def RiAttributeEnd():
    rlib.RiAttributeEnd()

def RiBound(bound ):
    rlib.RiBound(RtBound(*bound ))

def RiColor( *color ):
    rlib.RiColor( RtColor(*to_rfa(color) ))

def RiOpacity(*color):
    rlib.RiOpacity(RtColor(*to_rfa(color)))

def RiOption(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiOptionV(tok2c(name), n, toks, parms)

def RiReverseOrientation():
    rlib.RiReverseOrientation()

def RiTextureCoordinates(s1, t1, s2, t2, s3, t3, s4, t4):
    RiTextureCoordinates(to_rf(s1), to_rf(t1),
                         to_rf(s2), to_rf(t2),
                         to_rf(s3), to_rf(t3),
                         to_rf(s4), to_rf(t4))

def RiSides(sides):
    rlib.RiSides(ctypes.c_int(sides))

def RiDeclare(name, declaration):
    rlib.RiDeclar(tok2c(name), tok2c(declaration))

def RiLightSource(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    return rlib.RiLightSourceV(tok2c(name), n, toks, parms)

def RiAreaLightSource(name, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    return rlib.RiAreaLightSourceV(tok2c(name), n, toks, parms)

# TODO: Check RtBasis type is working here
# DL_INTERFACE RtVoid RiBasis(RtBasis ubasis, RtInt ustep,
# 						 RtBasis vbasis, RtInt vstep);
def RiBasis(ubasis, ustep, vbasis, vstep):
    rlib.RiBasis(RtBasis(*ubasis), ustep, RtBasis(*vbasis), vstep)

def RiPatch(ptype, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiPatchV(tok2c(ptype), n, toks, parms)

def RiPatchMesh(ptype, nu, uwrap, nv, vwrap, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiPatchMeshV(tok2c(ptype), to_ri(nu), tok2c(uwrap),
                      to_ri(nv), tok2c(vwrap), n, toks, parms)
                    
def RiPoints(npoints, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiPointsV(to_ri(npoints), n, toks, parms)

def RiCurves(ctype, ncurves, nvertices, wrap, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiCurvesV(tok2c(ctype), to_ri(ncurves), to_ri(nvertices), tok2c(wrap), n, toks, parms)


def RiNuCurves(ncurves, nvertices, order, knot, cmin, cmax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiNuCurvesV(to_ri(ncurves), to_ria(nvertices),
                     to_ria(order), to_rfa(knot),
                     to_rfa(cmin), to_rfa(cmax), n, toks, parms)

def RiNuPatch(nu, uorder, uknot, umin, umax,
              nv, vorder, vknot, vmin, vmax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiNuPatchV(to_ri(nu), to_ri(uorder), to_rfa(uknot), to_rf(umin), to_rf(umax),
                     to_ri(nv), to_ri(vorder), to_rfa(vknot), to_rf(vmin), to_rf(vmax), n, toks, parms)
              
def RiTrimCurve(nloops, ncurves, order, knot, cmin, cmax, n, u, v, w):
    rlib.RiTrimCurve(to_ri(nloops), to_ria(ncurves), to_ria(order), to_rfa(knot), to_rfa(cmin), to_rfa(cmax),
                     to_ria(n), to_rfa(u), to_rfa(v), to_rfa(w))

def RiSubdivisionMesh(scheme, nfaces, nvertices, vertices, ntags, tags, nargs, intargs, floatargs):
    (nn, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiSubdivisionMeshV(tok2c(scheme), to_ri(nfaces), to_ria(nvertices), to_ria(vertices),
                                        to_ri(ntags), to_rsa(tags), to_ria(nargs),
                                        to_ria(intargs), to_rfa(floatargs),
                                        nn, toks, parms)
def RiHierarchicalSubdivisionMesh(scheme, nfaces, nvertices, vertices, ntags,
                                  tags, nargs, intargs, floatargs, stringargs, *args, **kwargs):
    (nn, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiHierarchicalSubdivisionMeshV(tok2c(scheme), to_ri(nfaces), to_ria(nvertices), to_ria(vertices),
                                        to_ri(ntags), to_c(tags), to_ria(nargs),
                                        to_ria(intargs), to_rfa(floatargs), to_rsa(stringargs),
                                        nn, toks, parms)

def RiCone(height, radius, thetamax, *args, **kwargs):
    (nn, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiConeV(to_rf(height), to_rf(radius), to_rf(thetamax), nn, toks, parms)

def RiCylinder(height, radius, zmin, zmaxthetamax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiCylinderV(to_rf(height), to_rf(radius), to_rf(zmin), to_rf(zmax), to_rf(thetamax), n, toks, parms)

def RiDisk(height, radius, thetamax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiDiskV(to_rf(height), to_rf(radius), to_rf(theatmax), n, toks, parms)

def RiHyperboloid(p1, p2, thetamax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiHyperboloidV(RtPoint(*p1), RtPoint(*p2), to_rf(thetamax), n, toks, parms)

def RiParaboloid(rmax, zmin, zmax, thetamax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiParaboloidV(to_rf(rmax), to_rf(zmin), to_rf(zmax), to_rf(thetamax), n, toks, parms)

def RiSphere(radius,zmin,zmax, thetamax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiSphereV(to_rf(radius), to_rf(zmin), to_rf(zmax),
                   to_rf(thetamax), n, toks,parms)

def RiTorus(majorradius, minorradius, phimin, phimax, thetamax, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiTorusV(to_rf(majorradius), to_rf(minorradius), to_rf(phimin), to_rf(phimax), to_rf(thetamax), n, toks, parms)

def RiGeneralPolygon(nloops, nvertices, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiGeneralPolygonV(to_ri(nloops), to_ria(nvertices), n, toks, parms)

def RiBlobby(nleaf, nentry, entry, nfloat, floats, nstring, strings, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiBlobbyV(to_ri(nleaf), to_ri(nentry), to_ria(entry), to_ri(nfloat), to_rfa(floats), to_ri(nstring), to_rsa(strings), n, toks, parms)

def RiPointsGeneralPolygons(npolys, nloops, nvertices, vertices, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiPointsGeneralPolygonsV(to_ri(npolys), to_ria(nloops), to_ria(nvertices), to_ria(vertices), n, toks, parms)

def RiPointsPolygons(npolys, nvertices, vertices, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiPointsPolygonsV(to_ri(npolys), to_ria(nvertices), to_ria(vertices), n, toks, parms)


def RiPolygon(nvertices, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiPolygonV(to_ri(nvertices), n, toks, parms)

def RiColorSamples(n, nRGB, RGBn):
    rlib.RiColorSamples(to_ri(n), to_rfa(nRGB), to_rfa(RGBn))

def RiHider(htype, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiHiderV(tok2c(htype), n, toks, parms)

def RiDetail(*bound):
    rlib.RiDetail(RtBound(*bound))

def RiDetailRange(minvisible, lowertransition, uppertransition, maxvisible):
    rlib.RiDetailRange(to_rf(minvisible), to_rf(lowertransition), to_rf(uppertransition), to_rf(maxvisible))

def RiGeometricApproximation(gatype, value):
    rlib.RiGeometricApproximation(tok2c(gatype), to_rf(value))

def RiGeometry(gtype, *args, **kwargs):
    (n, toks, parms) = parse_extra_args(args, kwargs)
    rlib.RiGeometryV(tok2c(gtype), n, toks, parms)

def RiOrientation(orientation):
    rlib.RiOrientation(tok2c(orientation))

# TODO: implement RiProcedureal
# DL_INTERFACE RtVoid RiProcedural(RtPointer i_data, RtBound i_bound,
# 							  RtVoid (*i_Subdivfunc)(RtPointer, RtFloat),
# 							  RtVoid (*i_Freefunc)(RtPointer));
def RiRelativeDetail(relativedetail):
    rlib.RiRelativeDetail(to_rf(relativedetail))

# TODO: implement RiProcDelayedReadArchive
# DL_INTERFACE RtVoid RiProcDelayedReadArchive(RtPointer data,
# 										  RtFloat detail);

# TODO: implement RiProcDynamicLoad
# def RiProcDynamicLoad(RtPointer data, RtFloat detail):
#     rlib.RiProcDynamicLoad(RtPointer data, RtFloat detail)

# TODO: implement RiProcRunProgram
# def RiProcRunProgram(RtPointer data, RtFloat detail):
#     rlib.RiProcRunProgram(RtPointer data, RtFloat detail)

# /*
# 	WARNING: On windows, you should provide your own RiProcFree if you allocate
# 	the memory yourself. It is not safe to allocate in a module (exe or dll)
# 	and free in another.
# */
# def RiProcFree(RtPointer data):
#     rlib.RiProcFree(RtPointer data)

# TODO: Implement RiMake* functions
# DL_INTERFACE RtVoid RiMakeBump(
# 		const char *picturename, const char *texturename,
# 		RtToken swrap, RtToken twrap,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth, RtFloat twidth, ... );
# DL_INTERFACE RtVoid RiMakeBumpV(
# 		const char *picturename, const char *texturename,
# 		RtToken swrap, RtToken twrap,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth, RtFloat twidth,
# 		RtInt n, __RI_CONST RtToken tokens[],
# 		RtPointer parms[] );
# def RiMakeBump(picturename, texturename, swrap, twrap,filterfunc, swidth, twidth, *args, **kwargs):
#     (n, toks, parms) = parse_extra_args(args, kwargs)
#     rlib.RiMakeBumpV(tok2c(picturename), tok2c(texturename), tok2c(swrap), tok2c(twrap),

# DL_INTERFACE RtVoid RiMakeCubeFaceEnvironment(
# 		const char *px, const char *nx,
# 		const char *py, const char *ny,
# 		const char *pz, const char *nz,
# 		const char *texturename,
# 		RtFloat fov,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth,
# 		RtFloat twidth, ... );
# DL_INTERFACE RtVoid RiMakeCubeFaceEnvironmentV(
# 		const char *px, const char *nx,
# 		const char *py, const char *ny,
# 		const char *pz, const char *nz,
# 		const char *texturename,
# 		RtFloat fov,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth,
# 		RtFloat twidth,
# 		RtInt n,
# 		__RI_CONST RtToken tokens[],
# 		RtPointer parms[] );
# DL_INTERFACE RtVoid RiMakeLatLongEnvironment(
# 		const char *picturename,
# 		const char *texturename,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth,
# 		RtFloat twidth, ... );
# DL_INTERFACE RtVoid RiMakeLatLongEnvironmentV(
# 		const char *picturename,
# 		const char *texturename,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth,
# 		RtFloat twidth,
# 		RtInt n,
# 		__RI_CONST RtToken tokens[],
# 		RtPointer parms[] );
# DL_INTERFACE RtVoid RiMakeShadow(
# 		const char *picturename,
# 		const char *texturename, ... );
# DL_INTERFACE RtVoid RiMakeShadowV(
# 		const char *picturename,
# 		const char *texturename,
# 		RtInt n,
# 		__RI_CONST RtToken tokens[],
# 		RtPointer parms[] );
# DL_INTERFACE RtVoid RiMakeTexture(
# 		const char *picturename, const char *texturename,
# 		RtToken swrap, RtToken twrap,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth, RtFloat twidth, ... );
# DL_INTERFACE RtVoid RiMakeTextureV(
# 		const char *picturename, const char *texturename,
# 		RtToken swrap, RtToken twrap,
# 		RtFilterFunc filterfunc,
# 		RtFloat swidth, RtFloat twidth,
# 		RtInt n, __RI_CONST RtToken tokens[],
# 		RtPointer parms[] );

# DL_INTERFACE RtVoid RiMakeBrickMap(
# 		RtInt nptc, const char *const *ptcnames, const char *bkmname, ...);
# DL_INTERFACE RtVoid RiMakeBrickMapV(
# 		RtInt nptc, const char *const *ptcnames, const char *bkmname,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);

def RiBoxFilter(x, y, xwidth, ywidth):
    return float(rlib.RiBoxFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiTriangleFilter(x, y, xwidth, ywidth):
    return float(rlib.RiTriangleFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiCatmullRomFilter(x, y, xwidth, ywidth):
    return float(rlib.RiCatmullRomFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiSeparableCatmullRomFilter(x, y, xwidth, ywidth):
    return float(rlib.RiSeparableCatmullRomFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiBesselFilter(x, y, xwidth, ywidth):
    return float(rlib.RiBesselFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiGaussianFilter(x, y, xwidth, ywidth):
    return float(rlib.RiGaussianFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiSincFilter(x, y, xwidth, ywidth):
    return float(rlib.RiSincFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiBlackmanHarrisFilter(x, y, xwidth, ywidth):
    return float(rlib.RiBlackmanHarrisFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiMitchellFilter(x, y, xwidth, ywidth):
    return float(rlib.RiMitchellFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))

def RiDiskFilter(x, y, xwidth, ywidth):
    return float(rlib.RiDiskFilter(to_rf(x),to_rf(y),to_rf(xwidth), to_rf(ywidth)))


# TODO: Implement RiErrorHandler
# def RiErrorHandler( RtErrorHandler handler):
#     rlib.RiErrorHandler( RtErrorHandler handler)

def RiErrorIgnore(code, severity, msg):
    rlib.RiErrorIgnore(to_ri(code), to_ri(severity), tok2c(msg))

def RiErrorPrint(code, severity, msg):
    rlib.RiErrorPrint(to_ri(code), to_ri(severity), tok2c(msg))

def RiErrorAbort(code, severity, msg):
    rlib.RiErrorAbort(to_ri(code), to_ri(severity), tok2c(msg))

def RiSynchronize(tok):
    rlib.RiSynchronize(tok2c(tok))



# Error values
RIE_NOERROR     = 0

RIE_NOMEM       = 1
RIE_SYSTEM      = 2
RIE_NOFILE      = 3
RIE_BADFILE     = 4
RIE_VERSION     = 5
RIE_DISKFULL    = 6

RIE_INCAPABLE   = 11
RIE_UNIMPLEMENT = 12
RIE_LIMIT       = 13
RIE_BUG         = 14

RIE_NOTSTARTED  = 23
RIE_NESTING     = 24
RIE_NOTOPTIONS  = 25
RIE_NOTATTRIBS  = 26
RIE_NOTPRIMS    = 27
RIE_ILLSTATE    = 28
RIE_BADMOTION   = 29
RIE_BADSOLID    = 30

RIE_BADTOKEN    = 41
RIE_RANGE       = 42
RIE_CONSISTENCY = 43
RIE_BADHANDLE   = 44
RIE_NOSHADER    = 45
RIE_MISSINGDATA = 46
RIE_SYNTAX      = 47
RIE_TOKENREDECLARED = 48

RIE_MATH        = 61

RIE_BADATTRIB = 140
RIE_BADOPTION = 141
RIE_SPACEREDECLARED =142
RIE_NODISPLAY   = 143
RIE_ERRRERTOOBID  = 144
RIE_ERRBADSHADERPARAM = 145
RIE_ERRSHADERPARAMMISMATCH = 146
RIE_ERRBADARRAYACCESSINSHADER = 147

RIE_SHADER_PRINTF = 199

# Error severities
RIE_INFO        = 0
RIE_WARNING     = 1
RIE_ERROR       = 2
RIE_SEVERE      = 3
