from cadquery import *
from cqplugin import export_step
import unittest
import sys
import os

def makeCube(size, xycentered=True):
    if xycentered:
        return Workplane().rect(size, size).extrude(size).val()
    else:
        return Solid.makeBox(size, size, size)

def makeUnitCube(centered=True):
    return makeCube(1.0, centered)

def readFileAsString(fileName):
    f = open(fileName, "r")
    s = f.read()
    f.close()
    return s



