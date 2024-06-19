from cadquery import *
from cadquery.occ_impl.exporters.assembly import exportAssembly
from cqplugin.export_step import export_step_ap214
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
        
        #starting variables
        cwd = os.getcwd()
        tmpfile = "temp.step"
        tmp_path = os.path.join(cwd,tmpfile)

        #create something
        modell = makeUnitCube()
        assembly = Assembly()
        assembly.add(modell,loc=Location((0, 0, 0), (1, 0, 0), 180),name="Test text.")
        
        #save step to tmp file
        assembly.save(tmp_path)
        
        #read tmp file and comper result
        FileAsString = readFileAsString(tmp_path)
        
        #remove tmp file
        if os.path.exists(tmp_path): os.remove(tmp_path)
        
        self.assertRegex(FileAsString,"Test text.")

class TestSTEPap214export(unittest.TestCase):
    def test_step_ap214_export(self):
        cwd = os.getcwd()
        tmpfile = "temp.step"
        tmp_path = os.path.join(cwd,tmpfile)
        
        #create something
        modell = makeUnitCube()
        assembly = Assembly()
        assembly.add(modell,loc=Location((0, 0, 0), (1, 0, 0), 180))
        
        #save step to tmp file
        export_step_ap214(assembly,tmp_path)
        
        #read tmp file and comper result
        FileAsString = readFileAsString(tmp_path)
        
        #remove tmp file
        if os.path.exists(tmp_path): os.remove(tmp_path)
        
        self.assertRegex(FileAsString,"Test text.")
