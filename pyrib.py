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
import os

delight_path = os.getenv("DELIGHT")
delight_lib = delight_path + '/lib/lib3delight.dylib'
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

def tok2c(val):
    return ctypes.c_char_p(val.encode('utf8'))



def RiBegin(name=''):
    if len(name)==0:
        rlib.RiBegin(0)
    else:
        rlib.RiBegin(tok2c(cname))

def RiEnd():
    rlib.RiEnd()

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
# DL_INTERFACE RtVoid RiMotionBeginV(RtInt n, RtFloat times[]);

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
def RiObjectBegin():
    return rlib.RiObjectBegin()

# DL_INTERFACE RtObjectHandle RiObjectBeginV(
# 	RtInt n, __RI_CONST RtToken tokens[], RtPointer params[] );
def RiObjectBeginV(n, tokens, params):
    pass

def RiObjectEnd():
    rlib.RiObjectEnd()


def RiObjectInstance(handle):
    rlib.RiObjectInstance(handle)


# DL_INTERFACE RtVoid RiResource(RtToken handle, RtToken type, ...);
def RiResource(handle, rtype, *args):
    print("RiResource not supported yet!\n")
    pass

# DL_INTERFACE RtVoid RiResourceV(
# 		RtToken handle, RtToken type,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);
def RiResourceV(handle, rtype, n, tokens, params):
    print("RiResourceV not supported yet!\n")
    pass

def RiResourceBegin():
    rlib.RiResourceBegin()

def RiResourceEnd():
    rlib.RiResourceEnd()

def RiFormat(xres, yres, aspect):
    rlib.RiFormat(ctypes.c_int(xres), ctypes.c_int(yres), rf(aspect))
    
def RiFrameAspectRatio(aspect):
    rlib.RiFrameAspectRatio(rf(aspect))

def RiScreenWindow(left, right, bottom, top):
    rlib.RiScreenWindow(rf(left), rf(right), rf(bottom), rf(top))

def RiClipping(hither, yon):
    rlib.RiClipping(rf(hither), rf(yon))

def RiClippingPlane(x, y, z, nx, ny, nz ):
    rlib.RiClippingPlane(rf(x), rf(y), rf(z), rf(nx), rf(ny), rf(nz))

def RiCropWindow(xmin, xmax, ymin, ymax):
    rlib.RiCropWindow(rf(xmin), rf(xmax), rf(ymin), rf(ymax))
    
def RiDepthOfField(fstop, focallength, focaldistance):
    rlib.RiDepthOfField(rf(fstop), rf(focallength), rf(focaldistance))

def RiProjection(name, *args):
    toks = []
    parms = []
    for i in range(len(args)//2):
        toks.append(tok2c(args[i*2]))
        parms.append(0)
    
    rlib.RiProjectionV(tok2c(name), 0, 0,0)

# DL_INTERFACE RtVoid RiProjectionV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
def RiShutter(smin, max):
    rlib.RiShutter(rf(smin), rf(smax))

# def RiCamera(RtToken camera, ...):
#     rlib.RiCamera(RtToken camera, ...)

# DL_INTERFACE RtVoid RiCameraV(
# 		RtToken camera,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

def RiDisplay(name, dtype, mode, *args):
    toks = []
    parms = []
    for i in range(len(args)//2):
        toks.append(tok2c(args[i*2]))
        parms.append(0)
    rlib.RiDisplayV(tok2c(name), tok2c(dtype), tok2c(mode), ctypes.c_int(len(args)//2),0,0)
                  
		# const char *name, RtToken type, RtToken mode, ...);

# DL_INTERFACE RtVoid RiDisplayV(
# 		const char *name, RtToken type, RtToken mode,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiDisplayChannel(const char *channel, ...):
#     rlib.RiDisplayChannel(const char *channel, ...)

# DL_INTERFACE RtVoid RiDisplayChannelV(
# 		const char *channel,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);

def RiExposure(gain, gamma):
    rlib.RiExposure(rf(gain), rf(gamma))

# def RiImager(RtToken name, ...):
#     rlib.RiImager(RtToken name, ...)

# DL_INTERFACE RtVoid RiImagerV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiPixelFilter(RtFilterFunc filterfunc,
# 							   RtFloat xwidth, RtFloat ywidth);
def RiPixelSamples(xsamples, ysamples):
    rlib.RiPixelSamples(rf(xsamples), rf(ysamples))

def RiPixelVariance(variation):
    rlib.RiPixelVariance(rf(variation))

def RiQuantize(qtype, one, qmin, qmax, ampl):
    rlib.RiQuantize(tok2c(qtype), ctypes.c_int(one), ctypes.c_int(qmin), ctypes.c_int(qmax), rf(ampl));

# def RiConcatTransform(RtMatrix transform):
#     rlib.RiConcatTransform(RtMatrix transform)

def RiCoordinateSystem(space):
    rlib.RiCoordinateSystem(tok2c(space))

def RiScopedCoordinateSystem(space):
    rlib.RiScopedCoordinateSystem(tok2c(space))

def RiCoordSysTransform(space):
    rlib.RiCoordSysTransform(tok2c(space))

def RiIdentity():
    rlib.RiIdentity()

def RiPerspective(fov):
    rlib.RiPerspective(ctypes.c_cloat(fov))

def RiRotate(angle, dx, dy, dz):
    rlib.RiRotate(rf(angle), rf(dx), rf(dy), rf(dz))
def RiScale(dx, dy, dz):
    rlib.RiScale(rf(dx), rf(dy), rf(dz))
def RiTranslate(dx, dy, dz):
    rlib.RiTranslate(rf(dx), rf(dy), rf(dz))

def RiSkew(angle, dx1, dy1, dz1, dx2, dy2, dz2):
    rlib.RiSkew(rf(angle), rf(dx1), rf(dy1), rf(dz1),
              rf(dx2), rf(dy2), rf(dz2))

# def RiTransform(RtMatrix transform):
#     rlib.RiTransform(RtMatrix transform)

def RiTransformBegin(void):
    rlib.RiTransformBegin()

def RiTransformEnd(void):
    rlib.RiTransformEnd()

# DL_INTERFACE RtPoint* RiTransformPoints(RtToken fromspace,
# 									 RtToken tospace,
# 									 RtInt n, RtPoint points[]);

# def RiAtmosphere(RtToken name, ...):
#     rlib.RiAtmosphere(RtToken name, ...)

# DL_INTERFACE RtVoid RiAtmosphereV(
# 		RtToken name, RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiDeformation(RtToken name, ...):
#     rlib.RiDeformation(RtToken name, ...)

# DL_INTERFACE RtVoid RiDeformationV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiDisplacement(RtToken name, ...):
#     rlib.RiDisplacement(RtToken name, ...)

# DL_INTERFACE RtVoid RiDisplacementV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiExterior(RtToken name, ...):
#     rlib.RiExterior(RtToken name, ...)

# DL_INTERFACE RtVoid RiExteriorV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiIlluminate(RtLightHandle light, RtBoolean onoff):
#     rlib.RiIlluminate(RtLightHandle light, RtBoolean onoff)

# def RiInterior(RtToken name, ...):
#     rlib.RiInterior(RtToken name, ...)

# DL_INTERFACE RtVoid RiInteriorV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# def RiShader(RtToken name, RtToken handle, ...):
#     rlib.RiShader(RtToken name, RtToken handle, ...)

# DL_INTERFACE RtVoid RiShaderV( RtToken name, RtToken handle,
# 	RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

def RiMatte(onoff):
    rlib.RiMatte(ctypes.c_bool(onoff))

def RiMultiplyShadingRate(ratemultiplier):
    rlib.RiMultiplyShadingRate(rf(ratemultiplier))

def RiShadingRate(size):
    rlib.RiShadingRate(rf(size))

def RiShadingInterpolation(itype):
    rlib.RiShadingInterpolation(tok2c(itype))

# def RiSurface(RtToken name, ...):
#     rlib.RiSurface(RtToken name, ...)

# DL_INTERFACE RtVoid RiSurfaceV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# def RiArchiveRecord(RtToken type, __RI_CONST char *format, ...):
#     rlib.RiArchiveRecord(RtToken type, __RI_CONST char *format, ...)

# DL_INTERFACE RtVoid RiReadArchive(
# 		RtToken name, RtArchiveCallback callback, ...);
# DL_INTERFACE RtVoid RiReadArchiveV(
# 		RtToken name, RtArchiveCallback callback,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtArchiveHandle RiArchiveBegin(RtToken archivename, ...);
# DL_INTERFACE RtArchiveHandle RiArchiveBeginV(
# 		RtToken archivename,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);
def RiArchiveEnd():
    rlib.RiArchiveEnd()


# def RiIfBegin(RtToken expression, ...):
#     rlib.RiIfBegin(RtToken expression, ...)

# DL_INTERFACE RtVoid RiIfBeginV(
# 		RtToken expression,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);
# def RiElseIf(RtToken expression, ...):
#     rlib.RiElseIf(RtToken expression, ...)

# DL_INTERFACE RtVoid RiElseIfV(
# 		RtToken expression,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);
def RiElse():
    rlib.RiElse()

def RiIfEnd():
    rlib.RiIfEnd()

# def RiAttribute(RtToken name, ...):
#     rlib.RiAttribute(RtToken name, ...)

# DL_INTERFACE RtVoid RiAttributeV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
def RiAttributeBegin(void):
    rlib.RiAttributeBegin()

def RiAttributeEnd(void):
    rlib.RiAttributeEnd()

def RiBound(bound ):
    rlib.RiBound(RtBound(*bound ))

def RiColor( color ):
    rlib.RiColor( RtColor(*color) )

def RiOpacity(color):
    rlib.RiOpacity(RtColor(*color))

# def RiOption(RtToken name, ...):
#     rlib.RiOption(RtToken name, ...)

# DL_INTERFACE RtVoid RiOptionV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
def RiReverseOrientation():
    rlib.RiReverseOrientation()

def RiTextureCoordinates(s1, t1, s2, t2, s3, t3, s4, t4):
    RiTextureCoordinates(rf(s1), rf(t1),
                         rf(s2), rf(t2),
                         rf(s3), rf(t3),
                         rf(s4), rf(t4))

def RiSides(sides):
    rlib.RiSides(ctypes.c_int(sides))


def RiDeclare(name, declaration):
    rlib.RiDeclar(tok2c(name), tok2c(declaration))

# DL_INTERFACE RtLightHandle RiLightSource(RtToken name, ...);
# DL_INTERFACE RtLightHandle RiLightSourceV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtLightHandle RiAreaLightSource(RtToken name, ...);
# DL_INTERFACE RtLightHandle RiAreaLightSourceV(
# 		RtToken name,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtVoid RiBasis(RtBasis ubasis, RtInt ustep,
# 						 RtBasis vbasis, RtInt vstep);
# def RiPatch(RtToken type, ...):
#     rlib.RiPatch(RtToken type, ...)

# DL_INTERFACE RtVoid RiPatchMesh(RtToken type, RtInt nu, RtToken uwrap,
# 							 RtInt nv, RtToken vwrap, ...);
# DL_INTERFACE RtVoid RiPatchV(
# 		RtToken type,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiPatchMeshV(
# 		RtToken type,
# 		RtInt nu, RtToken uwrap,
# 		RtInt nv, RtToken vwrap,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiPoints(RtInt npoints, ...):
#     rlib.RiPoints(RtInt npoints, ...)

# DL_INTERFACE RtVoid RiPointsV(
# 		RtInt npoints,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiCurves(RtToken type, RtInt ncurves,
# 						  RtInt nvertices[], RtToken wrap, ...);
# DL_INTERFACE RtVoid RiCurvesV(
# 		RtToken type, RtInt ncurves,
# 		RtInt nvertices[], RtToken wrap,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiNuCurves(
# 		RtInt ncurves, RtInt nvertices[],
# 		RtInt order[], RtFloat knot[],
# 		RtFloat min[], RtFloat max[], ...);
# DL_INTERFACE RtVoid RiNuCurvesV(
# 		RtInt ncurves, RtInt nvertices[],
# 		RtInt order[], RtFloat knot[],
# 		RtFloat min[], RtFloat max[],
# 		RtInt n, __RI_CONST RtToken tokens[],
# 		RtPointer params[]);

# DL_INTERFACE RtVoid RiNuPatch(RtInt nu, RtInt uorder, RtFloat uknot[],
# 						   RtFloat umin, RtFloat umax,
# 						   RtInt nv, RtInt vorder, RtFloat vknot[],
# 						   RtFloat vmin, RtFloat vmax, ...);
# DL_INTERFACE RtVoid RiNuPatchV(RtInt nu, RtInt uorder, RtFloat uknot[],
# 							RtFloat umin, RtFloat umax,
# 							RtInt nv, RtInt vorder, RtFloat vknot[],
# 							RtFloat vmin, RtFloat vmax, RtInt n,
# 							__RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiTrimCurve(RtInt nloops, RtInt ncurves[],
# 							 RtInt order[], RtFloat knot[],
# 							 RtFloat min[], RtFloat max[], RtInt n[],
# 							 RtFloat u[], RtFloat v[], RtFloat w[]);

# DL_INTERFACE RtVoid RiSubdivisionMesh(
# 		RtToken scheme,
# 		RtInt nfaces, RtInt nvertices[], RtInt vertices[],
# 		RtInt ntags, __RI_CONST RtToken tags[], RtInt nargs[],
# 		RtInt intargs[], RtFloat floatargs[],
# 		...);
# DL_INTERFACE RtVoid RiSubdivisionMeshV(
# 		RtToken scheme,
# 		RtInt nfaces, RtInt nvertices[], RtInt vertices[],
# 		RtInt ntags, __RI_CONST RtToken tags[], RtInt nargs[],
# 		RtInt intargs[], RtFloat floatargs[],
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtVoid RiHierarchicalSubdivisionMesh(
# 		RtToken scheme,
# 		RtInt nfaces, RtInt nvertices[], RtInt vertices[],
# 		RtInt ntags, __RI_CONST RtToken tags[], RtInt nargs[],
# 		RtInt intargs[], RtFloat floatargs[], __RI_CONST RtToken stringargs[],
# 		...);
# DL_INTERFACE RtVoid RiHierarchicalSubdivisionMeshV(
# 		RtToken scheme,
# 		RtInt nfaces, RtInt nvertices[], RtInt vertices[],
# 		RtInt ntags, __RI_CONST RtToken tags[], RtInt nargs[],
# 		RtInt intargs[], RtFloat floatargs[], __RI_CONST RtToken stringargs[],
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer params[]);

# DL_INTERFACE RtVoid RiCone(RtFloat height, RtFloat radius,
# 						RtFloat thetamax, ...);
# DL_INTERFACE RtVoid RiConeV(
# 		RtFloat height, RtFloat radius,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiCylinder(RtFloat radius,RtFloat zmin,RtFloat zmax,
# 							RtFloat thetamax,...);
# DL_INTERFACE RtVoid RiCylinderV(
# 		RtFloat radius,RtFloat zmin,RtFloat zmax,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiDisk(RtFloat height, RtFloat radius,
# 						RtFloat thetamax, ...);
# DL_INTERFACE RtVoid RiDiskV(
# 		RtFloat height, RtFloat radius,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiHyperboloid(RtPoint point1, RtPoint point2,
# 							   RtFloat thetamax, ...);
# DL_INTERFACE RtVoid RiHyperboloidV(
# 		RtPoint point1, RtPoint point2,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiParaboloid(RtFloat rmax,RtFloat zmin,RtFloat zmax,
# 							  RtFloat thetamax,...);
# DL_INTERFACE RtVoid RiParaboloidV(
# 		RtFloat rmax,RtFloat zmin,RtFloat zmax,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);
def RiSphere(radius,zmin,zmax, thetamax, *args):
    rlib.RiSphereV(rf(radius), rf(zmin), rf(zmax),
                   rf(thetamax), 0, 0,0)

# DL_INTERFACE RtVoid RiSphereV(
# 		RtFloat radius,RtFloat zmin,RtFloat zmax,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiTorus(RtFloat majorradius, RtFloat minorradius,
# 						 RtFloat phimin, RtFloat phimax,
# 						 RtFloat thetamax, ...);
# DL_INTERFACE RtVoid RiTorusV(
# 		RtFloat majorradius, RtFloat minorradius,
# 		RtFloat phimin, RtFloat phimax,
# 		RtFloat thetamax, RtInt n,
# 		__RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtVoid RiGeneralPolygon(RtInt nloops,
# 								  RtInt nvertices[], ...);
# DL_INTERFACE RtVoid RiGeneralPolygonV(
# 		RtInt nloops, RtInt nvertices[],
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtVoid RiBlobby(
# 	RtInt nleaf,
# 	RtInt nentry, RtInt entry[],
# 	RtInt nfloat, RtFloat floats[],
# 	RtInt nstring, RtString strings[], ...);
# DL_INTERFACE RtVoid RiBlobbyV(
# 	RtInt nleaf,
# 	RtInt nentry, RtInt entry[],
# 	RtInt nfloat, RtFloat floats[],
# 	RtInt nstring, RtString strings[],
# 	RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtVoid RiPointsGeneralPolygons(RtInt npolys, RtInt nloops[],
# 										 RtInt nvertices[],
# 										 RtInt vertices[], ...);
# DL_INTERFACE RtVoid RiPointsGeneralPolygonsV(
# 		RtInt npolys, RtInt nloops[],
# 		RtInt nvertices[],
# 		RtInt vertices[],
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# DL_INTERFACE RtVoid RiPointsPolygons(RtInt npolys, RtInt nvertices[],
# 								  RtInt vertices[], ...);
# DL_INTERFACE RtVoid RiPointsPolygonsV(
# 		RtInt npolys, RtInt nvertices[],
# 		RtInt vertices[],
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
# def RiPolygon(RtInt nvertices, ...):
#     rlib.RiPolygon(RtInt nvertices, ...)

# DL_INTERFACE RtVoid RiPolygonV(
# 		RtInt nvertices,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

# DL_INTERFACE RtVoid RiColorSamples(RtInt n, RtFloat nRGB[],
# 								RtFloat RGBn[]);
# def RiHider(RtToken type, ...):
#     rlib.RiHider(RtToken type, ...)

# DL_INTERFACE RtVoid RiHiderV(
# 		RtToken type,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);

def RiDetail(*bound):
    rlib.RiDetail(RtBound(*bound))

# DL_INTERFACE RtVoid RiDetailRange(RtFloat minvisible,
# 							   RtFloat lowertransition,
# 							   RtFloat uppertransition,
# 							   RtFloat maxvisible);
def RiGeometricApproximation(gatype, value):
    rlib.RiGeometricApproximation(tok2c(gatype), rf(value))

# def RiGeometry(RtToken type, ...):
#     rlib.RiGeometry(RtToken type, ...)

# DL_INTERFACE RtVoid RiGeometryV(
# 		RtToken type,
# 		RtInt n, __RI_CONST RtToken tokens[], RtPointer parms[]);
def RiOrientation(orientation):
    rlib.RiOrientation(tok2c(orientation))

# DL_INTERFACE RtVoid RiProcedural(RtPointer i_data, RtBound i_bound,
# 							  RtVoid (*i_Subdivfunc)(RtPointer, RtFloat),
# 							  RtVoid (*i_Freefunc)(RtPointer));
def RiRelativeDetail(relativedetail):
    rlib.RiRelativeDetail(rf(relativedetail))

# DL_INTERFACE RtVoid RiProcDelayedReadArchive(RtPointer data,
# 										  RtFloat detail);
# def RiProcDynamicLoad(RtPointer data, RtFloat detail):
#     rlib.RiProcDynamicLoad(RtPointer data, RtFloat detail)

# def RiProcRunProgram(RtPointer data, RtFloat detail):
#     rlib.RiProcRunProgram(RtPointer data, RtFloat detail)

# /*
# 	WARNING: On windows, you should provide your own RiProcFree if you allocate
# 	the memory yourself. It is not safe to allocate in a module (exe or dll)
# 	and free in another.
# */
# def RiProcFree(RtPointer data):
#     rlib.RiProcFree(RtPointer data)


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

# DL_INTERFACE RtFloat RiBoxFilter(RtFloat x, RtFloat y,
# 							  RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiTriangleFilter(RtFloat x, RtFloat y,
# 								   RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiCatmullRomFilter(RtFloat x, RtFloat y,
# 									 RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiSeparableCatmullRomFilter(RtFloat x, RtFloat y,
# 									 RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiBesselFilter(RtFloat x, RtFloat y,
# 								 RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiGaussianFilter(RtFloat x, RtFloat y,
# 								   RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiSincFilter(RtFloat x, RtFloat y, RtFloat xwidth, RtFloat ywidth);
# DL_INTERFACE RtFloat RiBlackmanHarrisFilter(RtFloat i_x, RtFloat i_y, RtFloat i_filterXWidth, RtFloat i_filteryWidth);
# DL_INTERFACE RtFloat RiMitchellFilter(RtFloat i_x, RtFloat i_y, RtFloat i_filterXWidth, RtFloat i_filteryWidth);
# DL_INTERFACE RtFloat RiDiskFilter(RtFloat i_x, RtFloat i_y, RtFloat i_filterXWidth, RtFloat i_filteryWidth);


# def RiErrorHandler( RtErrorHandler handler):
#     rlib.RiErrorHandler( RtErrorHandler handler)

# def RiErrorIgnore(RtInt code, RtInt severity, __RI_CONST char *msg):
#     rlib.RiErrorIgnore(RtInt code, RtInt severity, __RI_CONST char *msg)

# def RiErrorPrint(RtInt code, RtInt severity, __RI_CONST char *msg):
#     rlib.RiErrorPrint(RtInt code, RtInt severity, __RI_CONST char *msg)

# def RiErrorAbort(RtInt code, RtInt severity, __RI_CONST char *msg):
#     rlib.RiErrorAbort(RtInt code, RtInt severity, __RI_CONST char *msg)

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
