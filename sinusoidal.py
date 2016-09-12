#!/usr/bin/env python
#
#----------------------------------------------------------------------
# sinusoidal.py --- Program to specify vertices, edges,
# and faces of a droplet on a sinusoidal surface with a single scale
# of roughness.
#
# Written by Michael S. Bell
# msb5390@gmail.com
# August 18, 2015
#----------------------------------------------------------------------

import glob
import os
import subprocess
import sys
#from optparse import OptionParser
import argparse
import math
from mike_math import step_function, sign
from dprint import dprint
import fcntl

#usage = "%prog [options]"
parser = argparse.ArgumentParser(description='Create a Surface Evolver script.')
parser.add_argument("-A", action="store", dest="amplitude", type=float,
                    metavar="<amp>", help="Amplitude")
parser.add_argument("-l", action="store", dest="lambda1", type=float,
                    metavar="<lambda>", help="lambda1")
parser.add_argument("-n", action="store", dest="nGrooves", type=int,
                  metavar="<nGrooves>", help="number of grooves")
parser.add_argument("-M", action="store", dest="mode",
                    choices=["C","W"],
                  metavar="<wetting mode>", help="wetting mode")
parser.add_argument("-b", action="store", dest="Bo", type=float,
                  metavar="<bond number>", help="bond number")
parser.add_argument("-e", action="store", dest="thetaE", type=float,
                  metavar="<thetaE>", help="thetaE in degrees")
args = parser.parse_args()


# Set to True to print out helpful debugging messages.
# Set to False to turn them off.
debug = True

#-----------------------------------------------------------------
# Droplet modeling parameters:
#------------------------------------------------------------------
# Today's date:
date = "August 18, 2015"

gravQuant = "grav"
pi = math.pi
# Surface- (and drop-) specific parameters:
R0 = 1.0
V = math.pi*pow(R0,2)

if args.thetaE == None: thetaE = 115.0
else: thetaE = args.thetaE

if args.amplitude == None: amplitude = 0.5
else: amplitude = args.amplitude

if args.lambda1 == None: lambda1 = 1.0
else: lambda1 = args.lambda1

if args.nGrooves == None: nGrooves = 3
else: nGrooves = args.nGrooves

if args.mode == None: mode = "Wenzel"
else:
    if args.mode == "C":
        mode = "Cassie"
    else:
        mode = "Wenzel"

if args.Bo == None: Bo = 0.0
else: Bo = args.Bo



#------------------------------------------------------------------
# Some useful calculations.
#------------------------------------------------------------------
# Even or odd n: store +1 for even, store -1 for odd.
evenOdd = 1
if nGrooves % 2 == 1:
    evenOdd = -1
# Initial x-value of the contact line.
x0 = (2*nGrooves + 1)*(lambda1/4)
#if nGrooves == 0:
#  x0 = 1
# Compute the height of the drop above the top of the surface
# required to make a rectangular-top drop with the correct
# volume.
## First, define the additional height variable.
hAdd = 0
# Now, give it a value, but throw an exception if it is negative.
try:
    if mode == "Cassie":
        hAdd = (pi - amplitude*lambda1*((2*nGrooves + 1)/2 - (nGrooves + 1)/pi))/((2*nGrooves + 1)/2*lambda1)
    else:
        hAdd = (pi - amplitude*lambda1*((2*nGrooves + 1)/2 - 1/pi))/((2*nGrooves + 1)/2*lambda1)
    # If hAdd is negative, throw an exception.
    if hAdd < 0:
        raise Exception('nGroovesException')
# Catch any exception that might be thrown.
except Exception as inst:
    print(type(inst))
    print(inst + ": trying to make drop span too many grooves for the chosen drop size.")
    sys.exit("Fatal error occurred with drop specification. Closing.")

# Define a string containing the expression for the function
# defining the boundary.
surfaceExpression = "{eOdd}*amp*cos(2*pi*p1/lambda)+amp".format(eOdd=evenOdd)
surfaceExpression1 = "{eOdd}*amp*cos(2*pi*x/lambda)+amp".format(eOdd=evenOdd)
dprint(surfaceExpression, debug)


#------------------------------------------------------------------
# Program the constraints.
#------------------------------------------------------------------
# Some points must lie on or above the surface.
surfaceTopConstraint = """\
constraint surfaceTop nonnegative
formula: y - ({surfExpression})

""".format(surfExpression=surfaceExpression1)
dprint(surfaceTopConstraint, debug)

# Package all of the constraints together into a single constraints string.
constraints = surfaceTopConstraint
constraintList = "constraint surfaceTop"

#------------------------------------------------------------------
# Program the boundaries.
#------------------------------------------------------------------
# Boundary 1 will be for contact points on a downward slope.
boundary1 = """
boundary 1 parameters 1
x1: p1
x2: {surfExpression}
energy:
e1: -Bo*amp/(16*pi^2)*(4*{eOdd}*sin(tilt_angle*pi/180.0)*((2*pi^2*p1^2-lambda^2)\
*cos(2*pi*p1/lambda)+lambda*(lambda-2*pi*p1*sin(2*pi*p1/lambda)))\
+amp*pi*cos(tilt_angle*pi/180.0)*(4*(2+({eOdd})^2)*pi*p1\
+8*{eOdd}*lambda*sin(2*pi*p1/lambda)\
+({eOdd})^2*lambda*sin(4*pi*p1/lambda)))\
+SURFT*(lambda/(2*pi))*incompleteEllipticE(2*pi*x/lambda, -(2*pi*amp/lambda)^2)
content:
c1: (-1)*{eOdd}*amp*lambda/(2*pi)*sin(2*pi*p1/lambda)-amp*p1

""".format(surfExpression=surfaceExpression,eOdd=evenOdd)

# Boundary 2 will be for contact points on an upward slope.
boundary2 = """
boundary 2 parameters 1
x1: p1
x2: {surfExpression}
energy:
e1: -(-Bo*amp/(16*pi^2)*(4*{eOdd}*sin(tilt_angle*pi/180.0)*((2*pi^2*p1^2-lambda^2)\
*cos(2*pi*p1/lambda)+lambda*(lambda-2*pi*p1*sin(2*pi*p1/lambda)))\
+amp*pi*cos(tilt_angle*pi/180.0)*(4*(2+({eOdd})^2)*pi*p1\
+8*{eOdd}*lambda*sin(2*pi*p1/lambda)\
+({eOdd})^2*lambda*sin(4*pi*p1/lambda)))\
+SURFT*(lambda/(2*pi))*incompleteEllipticE(2*pi*x/lambda, -(2*pi*amp/lambda)^2))
content:
c1: -1*{eOdd}*amp*lambda/(2*pi)*sin(2*pi*p1/lambda)-amp*p1

""".format(surfExpression=surfaceExpression,eOdd=evenOdd)

# Boundary 3 will be for display only.
boundary3 = """
boundary 3 parameters 1
x1: p1
x2: {surfExpression}

""".format(surfExpression=surfaceExpression)

boundaries = boundary1 + boundary2 + boundary3
#------------------------------------------------------------------
# Define the gravitational potential to allow for tilting.
#------------------------------------------------------------------
potentialQuantity = """
quantity {grav} energy method edge_vector_integral global
vector_integrand:
q1: -y^2/2*Bo*cos(tilt_angle*pi/180.0)
q2: -x^2/2*Bo*sin(tilt_angle*pi/180.0)

""".format(grav=gravQuant)


# Define some strings that will be printed directly to the datafile
# (once they are filled):
vertices = ""
edges = ""
faces = ""
bodies = ""
#----------------------------------------------------------
# Fill the vertices, edges, and faces lists.
#----------------------------------------------------------

# First we will focus on the vertices on the boundary
# on the right half of the drop.

iMin = 0        # Where to start placing vertices on boundary.
vertCounter = 1 # Keep a running total of vertices.
edgeCounter = 1 # Keep a running total of edges.
faceCounter = 1 # Keep a running total of the faces.
boundaryVertices = [] # This will store the boundary vertices for reuse on left.
topVertices = []
botVertices = []
bInit = 1       # What is first boundary?
inc = 1         # A ticker to keep track of which boundary is next.
# If the drop is in Wenzel, first vertex is at far right end of drop.
if mode == "Wenzel":
    iMin = nGrooves

# Keep track of consecutive vertices. When it hits 2, we make an edge
# to connect them.
consecVerts = 1
joinLastToFirst = True

# If the drop spans an odd number of grooves, change the defaults we
# set earlier.
if evenOdd == -1:
    bInit = 2
    consecVerts = 0
    inc = -1
    joinLastToFirst = False    

if mode == "Wenzel":
    bInit = 1

boundary = bInit

# Loop over vertices on boundary on right half of drop.
for i in range(0, nGrooves+1):
    # Determine the current x-coordinate
    x = (2*i+1)*lambda1/4
    botVertices.append([x])
    if i == nGrooves:
        topVertices.append([x])
    if ((i == iMin) or (i > iMin)):
        # Add this boundary location to list for re-use on left side.
        boundaryVertices.append([x, boundary])
        # Make a string to write to datafile for current vertex.
        currentVertString = "{i}     {x} boundary {b}".format(i=vertCounter,
                                                              x=x,
                                                              b=boundary)
        # Append the vertex string to the master list of vertices.
        vertices += currentVertString + "\n"

        # If we have made 2 "consecutive" vertices, join them with an edge.
        if consecVerts == 2:
            # Write the edge string.
            currentEdgeString = "{i}     {v1} {v2} density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                          v1=vertCounter-1,
                                                                                          v2=vertCounter,
                                                                                          constraints=constraintList,
                                                                                          grav=gravQuant)
            # Add current edge to master edge list.
            edges += currentEdgeString + "\n"

            # We also need to define a face for this edge since it is disjoint.
            currentFaceString = "{i}     {edge}     {constraints}".format(i=faceCounter,
                                                                          edge=edgeCounter,
                                                                          constraints=constraintList)
            faces += currentFaceString + "\n"
            faceCounter += 1
        
            # Increment the edge counter.
            edgeCounter += 1
            # Restart the consecutive vertex counter (since we are on boundaries).
            consecVerts = 0
        boundary = boundary + inc
        inc = inc*(-1)
        vertCounter = vertCounter + 1
        # Increment the number of consecutive vertices.
        consecVerts += 1
        
# Make a list of y-coordinates for vertices on the side-walls
# of the drop.
#nVert = int((hAdd + amplitude)/(lambda1/2))
#nVert = int((hAdd + amplitude)/(lambda1))
nVert = 2
topFaceEdge1 = edgeCounter
rightVert = vertCounter - 1
vertVertices = []
for i in range(1, nVert):
    # Keep a list of the vertical vertices for use on left side.
    y = amplitude + i*(hAdd + amplitude)/nVert
    vertVertices.append([x0, y])
    currentVertString = "{i}     {x} {y}     {constraints}".format(i=vertCounter,
                                                                   x=x0,
                                                                   y=y,
                                                                   constraints=constraintList)
    vertices += currentVertString + "\n"
    currentEdgeString = "{i}     {v1} {v2}     density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                      v1=vertCounter-1,
                                                                                      v2=vertCounter,
                                                                                      constraints=constraintList,
                                                                                      grav=gravQuant)
    edges += currentEdgeString + "\n"
    vertCounter += 1
    edgeCounter += 1

# Now loop over the x-coordinates on the boundary to make vertices
# across the top of the drop.
for vert in reversed(topVertices):
    y = 2*amplitude + hAdd
    x = vert[0]
    currentVertString = "{i}     {x} {y}     {constraints}".format(i=vertCounter,
                                                                   x=x, y=y,
                                                                   constraints=constraintList)
    vertices += currentVertString + "\n"
    

    currentEdgeString = "{i}     {v1} {v2}     density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                      v1=vertCounter-1,
                                                                                      v2=vertCounter,
                                                                                      constraints=constraintList,
                                                                                      grav=gravQuant)
    edges += currentEdgeString + "\n"
    vertCounter += 1
    edgeCounter += 1

# Now continue across top of drop on left side.
for vert in topVertices:
    y = 2*amplitude + hAdd
    x = -vert[0]
    currentVertString = "{i}     {x} {y}     {constraints}".format(i=vertCounter,
                                                                   x=x, y=y,
                                                                   constraints=constraintList)
    vertices += currentVertString + "\n"

    currentEdgeString = "{i}     {v1} {v2}     density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                      v1=vertCounter-1,
                                                                                      v2=vertCounter,
                                                                                      constraints=constraintList,
                                                                                      grav=gravQuant)
    edges += currentEdgeString + "\n"
    vertCounter += 1
    edgeCounter += 1

for vert in reversed(vertVertices):
    x = -vert[0]
    y = vert[1]
    currentVertString = "{i}     {x} {y}     {constraints}".format(i=vertCounter,
                                                                   x=x, y=y,
                                                                   constraints=constraintList)
    vertices += currentVertString + "\n"

    currentEdgeString = "{i}     {v1} {v2}     density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                       v1=vertCounter-1,
                                                                                       v2=vertCounter,
                                                                                       constraints=constraintList,
                                                                                       grav=gravQuant)
    edges += currentEdgeString + "\n"
    vertCounter += 1
    edgeCounter += 1

leftVert = vertCounter - 1
# Finally, make the vertices on the left half of the base of the drop.
boundary = 2
inc = -1
consecVerts = -1
leftVert = 0
for vert in reversed(boundaryVertices):
    x = -vert[0]
    currentVertString = "{i}     {x} boundary {b}".format(i=vertCounter,
                                                          x=x,
                                                          b=boundary)
    # Append the vertex string to the master list of vertices.
    vertices += currentVertString + "\n"

    if consecVerts == -1:
        currentEdgeString = "{i}     {v1} {v2} density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                      v1=vertCounter-1,
                                                                                      v2=vertCounter,
                                                                                      constraints=constraintList,
                                                                                      grav=gravQuant)
        edges += currentEdgeString + "\n"
        edgeCounter += 1

        # Now make a face corresponding to all the top vertices.
        topFaceEdgeF = edgeCounter-1
        topFaceEdges = ""
        leftVert = vertCounter
        for i in range(topFaceEdge1, topFaceEdgeF+1):
            topFaceEdges += "{i} ".format(i=i)
        currentFaceString = "{i}     {edges}     {constraints}".format(i=faceCounter,
                                                                      edges=topFaceEdges,
                                                                      constraints=constraintList)
        faces += currentFaceString + "\n"
        faceCounter += 1
    
    # Increment the number of consecutive vertices.
    consecVerts += 1
    # If we have made 2 "consecutive" vertices, join them with an edge.
    if consecVerts == 2:
        # Write the edge string.
        currentEdgeString = "{i}     {v1} {v2} density 1 {constraints} {grav}".format(i=edgeCounter,
                                                                                      v1=vertCounter-1,
                                                                                      v2=vertCounter,
                                                                                      constraints=constraintList,
                                                                                      grav=gravQuant)
        # Add current edge to master edge list.
        edges += currentEdgeString + "\n"

        # We also need to define a face for this edge since it is disjoint.
        currentFaceString = "{i}     {edge}     {constraints}".format(i=faceCounter,
                                                                      edge=edgeCounter,
                                                                      constraints=constraintList)
        faces += currentFaceString + "\n"
        faceCounter += 1
        
        # Increment the edge counter.
        edgeCounter += 1
        # Restart the consecutive vertex counter (since we are on boundaries).
        consecVerts = 0
    boundary = boundary + inc
    inc = inc*(-1)
    vertCounter = vertCounter + 1

if consecVerts == 1:
    currentEdgeString = "{i}     {v1} {v2}     {constraints} {grav}".format(i=edgeCounter,
                                                                             v1=vertCounter-1,
                                                                             v2 = 1,
                                                                             constraints=constraintList,
                                                                             grav=gravQuant)
    edges += currentEdgeString + "\n"
    currentFaceString = "{i}     {edge}     {constraints}".format(i=faceCounter,
                                                                  edge=edgeCounter,
                                                                  constraints=constraintList)
    faces += currentFaceString + "\n"
    faceCounter+=1
    edgeCounter += 1


faceList = ""
for i in range(1,faceCounter):
    faceList += "{i} ".format(i=i)
bodies = "1     {faces} volume {V} density 1".format(faces=faceList,
                                                         V=V)



# Now, let's make the display surface.
vertices += """
// Display vertices:
"""

edges += """
// Display edges:
"""

t0 = -(nGrooves+1)*lambda1
displayEdge = 0
if evenOdd == -1:
    t0 += 0.5*lambda1

tf = -t0
nPts = int((tf - t0)/(lambda1/2))
for i in range(0,nPts+1):
    t = t0 + i*(lambda1/2)
    currentVertString = "{i}     {x} boundary 3 fixed".format(i=vertCounter,
                                                              x=t)
    vertices += currentVertString + "\n"
    if displayEdge != 0:
        currentEdgeString = "{i}     {v1} {v2} density 0 boundary 3 fixed noncontent".format(i=edgeCounter,
                                                                                             v1=vertCounter-1,
                                                                                             v2=vertCounter)
        edges += currentEdgeString + "\n"
        edgeCounter += 1
    vertCounter += 1
    displayEdge += 1


#------------------------------------------------------------------
# Comment block parameters (including program description):
#------------------------------------------------------------------
newfilename = "dropSinusoidalTest.fe"
# Description of the .fe file.
description = """// \
This program creates a SurfaceEvolver data file containing the vertices,
// edges and facets of a drop on a sinusoidal surface with specified
// dimensions.
// """
datafilePath = "/gpfs/home/azs5619/work/scripts/evolverset/"
myName = "Michael S. Bell"
myEmail = "msb5390@gmail.com"

#------------------------------------------------------------------
# Define parameters to be used by SurfaceEvolver
#------------------------------------------------------------------
PARAMETERS = """\
PARAMETER thetaE = {thetaE}
PARAMETER mode = \"{mode}\"
PARAMETER amp = {amp}
PARAMETER lambda = {lambda1}
PARAMETER tilt_angle = 0
PARAMETER Bo = {bond}""".format(thetaE=thetaE, mode=mode,
                           amp=amplitude, lambda1=lambda1, bond=Bo)
dprint(PARAMETERS, debug)


#------------------------------------------------------------------
# Define constants for the surface.
#------------------------------------------------------------------
definitions = """\
#define SURFT -cos({thetaE}*pi/180.0)
""".format(thetaE=thetaE)
dprint(definitions, debug)

#------------------------------------------------------------------
# Now construct the string that constitutes the SurfaceEvolver script.
#------------------------------------------------------------------
commentBlock = """//--------------------------------------------------------\
-------------------
// {filename}
{description} 
//
// {name}
// {email}
// {date}
//--------------------------------------------------------\
------------------""".format(filename=newfilename, description=description, name=myName, email=myEmail, date=date)

scriptText = """{commentBlock}

STRING
SPACE_DIMENSION 2

{PARAMETERS}

{definitions}
gravity_constant 0

{constraints}

{boundaries}

vertices:
{vertices}

edges:
{edges}

faces:
{faces}

bodies:
{bodies}

read
read "zebra.cmd"
//{{refine edges where on_boundary 3;}} 5
zebra
set edges where on_boundary 3 color black
window_aspect_ratio:=1

// Some commands to help with the evolution process.


// refine and update the zebra pattern
rz := {{ refine edges where not on_boundary 3; zebra;
	set edges where on_boundary 3 color black;
      }}

//error calling
call:={{printf "%3.9s \\n", "TILTANGLE";
  printf "%3.2f \\n", tilt_angle;
  printf "%3.9s \\n", "AMPLITUDE";
  printf "%3.2f \\n", amp;
  printf "%3.9s \\n", "BOND_NUMBER";
  printf "%3.2f \\n", bo;
  printf "%3.9s \\n", "LAMBDA";
  printf "%3.2f \\n", lambda;
  printf "%3.9s \\n", "THETA_E";
  printf "%3.2f \\n", thetaE;}}

//Function to get left vertex x-coordinate value of drop
left_vx :={{printf "%3.2f \\n", vertex[{lv}].x;}}

//Function to get right vertex x-coordinate value of drop
right_vx :={{printf "%3.2f \\n", vertex[{rv}].x;}}

//Function to get left vertex y-coordinate value of drop
left_vy :={{printf "%3.2f \\n", vertex[{lv}].y;}}

//Function to get right vertex y-coordinate value of drop
right_vy :={{printf "%3.2f \\n", vertex[{rv}].y;}}

// Function to return the apparent CA at right contact line.
car_app := {{if edge[{re}].x/edge[{re}].y > 0 then printf "%3.2f \\n", 180.0-atan(edge[{re}].y/edge[{re}].x)*180.0/pi
      else printf "%3.2f \\n", -atan(edge[{re}].y/edge[{re}].x)*180.0/pi;}}


// Function to return the actual CA at the right contact line.
car := {{if edge[{re}].x/edge[{re}].y > 0 then printf "%3.2f \\n", -atan({eOdd}*2*pi*amp/lambda*sin(2*pi*vertex[{rv}].x/lambda))*180/pi+180.0-atan(edge[{re}].y/edge[{re}].x)*180.0/pi
  else printf "%3.2f \\n", -atan({eOdd}*2*pi*amp/lambda*sin(2*pi*vertex[{rv}].x/lambda))*180/pi-atan(edge[{re}].y/edge[{re}].x)*180.0/pi;}}


// Function to return the apparent CA at the left contact line.
cal_app := {{foreach vertex[{lv}].edge do myEdge:=id;
  set edge[myEdge] color red;
  //print myEdge;
  if edge[myEdge].x/edge[myEdge].y > 0 then printf "%3.2f \\n", atan(edge[myEdge].y/edge[myEdge].x)*180.0/pi
        else printf "%3.2f \\n", 180 + atan(edge[myEdge].y/edge[myEdge].x)*180.0/pi;
  set edge[myEdge] color white;}}



// Function to return the actual CA at the left contact line.
cal := {{foreach vertex[{lv}].edge do myEdge:=id;
  set edge[myEdge] color red;
  //print myEdge;
  if edge[myEdge].x/edge[myEdge].y > 0 then printf "%3.2f \\n", atan({eOdd}*2*pi*amp/lambda*sin(2*pi*vertex[{lv}].x/lambda))*180/pi+atan(edge[myEdge].y/edge[myEdge].x)*180.0/pi
  else printf "%3.2f \\n", atan({eOdd}*2*pi*amp/lambda*sin(2*pi*vertex[{lv}].x/lambda))*180/pi+ 180 + atan(edge[myEdge].y/edge[myEdge].x)*180.0/pi;
  set edge[myEdge] color white;}}

go1 := {{g 50;
  rz;
  g 400;
  print "Now set Bo and tilt_angle!";}}

go2 := {{g 700;
  rz;
  g 13000;
  print "car = "; car;
  print "cal = "; cal;
  print "xA = "; print vertex[1].x;
  print "xR = "; print vertex[7].x;
  print "car_app = "; car_app;
  print "cal_app = "; cal_app;
  }}

fun := {{print "starting point of fun";
  Bo := 1.0; {{g10; V 100;}} 20; rz 2;      
  {{V 100; g 10;}} 25;
  print "finished fun";
  }}

function real average_finder(){{num_total := 0; 
  number_of_edges := 0;
  foreach edge where on_constraint surfaceTop do{{
    lengther:= edge[id].length;   
    printf "Edge Number  %d   Edge Length   %f\n", id, lengther;
    label:= edge[id].length;
    num_total:= num_total + label;
    number_of_edges := number_of_edges + 1;     
    
  }};
  final_answer := num_total/number_of_edges;
  return final_answer;

}}

function real standard_dev(){{
  stand_dev := 0;
  number_of_edges:=0;
  sumnum:=0;
  mean:= average_finder();
  foreach edge where on_constraint surfaceTop do{{
    stand_dev := stand_dev + sqr(edge[id].length - mean);
    number_of_edges := number_of_edges + 1;
  }};
  sumnum := 1 / (number_of_edges-1);
  stand_dev := sqrt(stand_dev/number_of_edges);
  return stand_dev;
}} 

stdev := {{printf "%f \n", standard_dev();}}
    

averageLength:={{
  summer:= average_finder();
  printf "%f \n", summer;
}}

function real hit_or_miss(real number){{
  functionvalue := number;
  op := 0;
  foreach vertex where on_constraint surfaceTop and value_of_constraint surfaceTop < functionvalue do{{
    op:= op + 1;
  }};
  return op;
}} 


// This "am_I_close" function checks if you are close to the sinusoidal surface;
// It will print True or "0" if you are too close to the surface;
// It will print False or "1" if you are fine and far away from the surface;

am_I_close:={{
  averager := average_finder();
  boolboolT := "True";
  boolboolF := "False";
  functional_output := 2*amp - (amp*cos(2*pi*averager/(2*lambda))+amp);
  whatsyournumber:= hit_or_miss(functional_output);
  if op = 0 then printf "%s \n", boolboolF else printf "%s \n", boolboolT;

  //printf "%f \n", whatsyournumber;
  //printf "%f \n", functional_output
}}

even_out:={{
  mean:= average_finder();
  foreach edge where on_constraint surfaceTop do{{
    while edge[id].length/mean > 1.25 do {V;};
    while edge[id].length/mean < 0.75 do {V;};
  }}
}} 


//End of Input
\
""".format(commentBlock=commentBlock, PARAMETERS=PARAMETERS,
           definitions=definitions,
           constraints=constraints+potentialQuantity,
           boundaries=boundaries,
           vertices=vertices, edges=edges,
           faces=faces, bodies=bodies,
           re=topFaceEdge1,
           rv=rightVert,
           le=topFaceEdgeF,
           lv=leftVert, 
	         eOdd=evenOdd)

fullfilepath = datafilePath + newfilename
feFile = open(fullfilepath, 'w') # open the file for writing (overwrite existing)
print(scriptText)
feFile.write(scriptText)
feFile.close()
