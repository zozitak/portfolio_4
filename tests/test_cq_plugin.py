from cadquery import *
from cqplugin import export_step
import unittest
import sys
import os

def makeUnitCube():
    return Workplane().rect(1.0, 1.0).extrude(1.0).val()

def readFileAsString(fileName):
    f = open(fileName, "r")
    s = f.read()
    f.close()
    return s

class TestSTEPexport(unittest.TestCase):
    def test_step_export(self):
        cwd = os.getcwd()
        tmpfile = "temp.step"
        #create something
        modell = makeUnitCube()
        assembly = Assembly()
        assembly.add(modell,loc=Location((0, 0, 0), (1, 0, 0), 180),name="Test text.")
        #save step to tmp file
        assembly.save(os.path.join(cwd,tmpfile))
        #read tmp file and comper result
        FileAsString = readFileAsString(os.path.join(cwd,tmpfile))
        #remove tmp file
        os.remove(os.path.join(cwd,tmpfile))
        self.assertRegex(FileAsString,"Test text.")

class TestSTEPap214export(unittest.TestCase):
    def test_step_ap214_export(self):
        raise NotImplemented
