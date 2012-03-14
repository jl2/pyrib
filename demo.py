#!/usr/bin/env python3

# demo.py

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


# A small demonstration of ri.py and some of its features

import math
import sys

import ri

def main(args):

    # Calls are very similar to the C API
    # RiBegin() -> ri.Begin()
    ri.Begin()

    ri.Display("demo.tif","file","rgba")
    ri.Format(800, 600,  1.0)

    ri.Projection("perspective")

    ri.Attribute("light", "string shadow", "on")

    # Variable arugment functions with token/argument pairs can be called using
    # two methods, both of which are shown here:
    ri.LightSource("distantlight",
                   # Technique 1, similar to C:
                   "from", [40.0,80.0,-40.0],
                   # Technique 2, using Python's keyword arguments
                   to=[0.0,0.0,0.0])

    # Three things to keep in mind:
    #    Some parameters, like "from" are Python keywords, so cannot
    #       be passed using keyword arguments
    #    Due to a Python limitation, all parameters using the keyword
    #       style must come after all of the parameters using the other style
    #    In many cases the ".0" in "0.0", "1.0", ... are important for autodetecting the C
    #       data type that is used
    #       "1" gets translated as the C type int while "1.0" is translated to an RtFloat
    #       Some situations are handled correctly without the ".0", but some cases still
    #       require data of the correc type be passed in.
    

    # Multiple word paramters must use the first parameter passing technique
    ri.Imager("background", "color background", [0.0, 0.25, 0.25])


    ri.Translate(0, 0,16)

    ri.WorldBegin()
        
    ri.AttributeBegin()
    ri.Color(1,0,0)
    ri.Surface("checker", blackcolor=[0.0,0.0,1.0])
    ri.Sphere(8, -8, 8, 360)
    ri.AttributeEnd()

    ri.AttributeBegin()
    ri.Color(0,1,0)
    ri.Surface("plastic")

    # Three ways of passing P, N, ... data:

    
    ri.Polygon(4, P =
               # "raw" values
               [10.0, -8.0, -10.0,
                # Grouped into tuples
                (10.0, -8.0, 10.0),
                # Grouped into lists
                [-10.0, -8.0, 10.0],
                -10.0, -8.0, -10.0,
                ])

    ri.AttributeEnd()
        
    ri.WorldEnd()

    ri.End();

if __name__=='__main__':
    main(sys.argv[1:])
