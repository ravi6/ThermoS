# Author: Ravi Saripalli
# 3rd Oct 2025
#Post process Heater Data
# and generate a report at a
# specific time stamp

import numpy as np
import sys

fname = "work/plant_res.csv"
with open(fname, 'r') as f:
   hdr = f.readline().strip().replace("\"","").split(',')

tCol = hdr.index ('time')
data = np.genfromtxt(fname, delimiter=',', skip_header=1)

def getVar (name, ts):
  idxVar = hdr.index (name.replace("\"","'"))
  var = np.interp (ts, data[:,tCol], data[:,idxVar])
  return (var) 

def report (uop, ts):
  inlet_m_flow = getVar (uop + ".inlet.m_flow", ts)
  Tf = getVar (uop + ".Tf", ts) - 273
  Tw = getVar (uop + ".Tw", ts) - 273
  Q_ew = getVar (uop + ".Q_ew", ts) / 1000
  Q_wf = getVar (uop + ".Q_wf", ts) / 1000
 # Cp   = getVar (uop + ".Cp", ts)
  dp = getVar(uop + ".inlet.p", ts) - 1e5 ;
  print ("Pressure Drop =", dp)

  print("=============<< Heater Summary >>=============")
  print(f"Tag:  {uop}  \t  Sample Time (s): {ts}")
  print("Inlet Conditions:")
  print(f"Flow rate (kg/s) = {inlet_m_flow:.8f}")
  print("Exit Conditions:")
  print(f"Fluid Temp (C) = {Tf:.1f}") 
  print(f"Wall Temp (C) = {Tw:.1f}")
  print(f"Q_ew (kW) = {Q_ew:.1f}")
  print(f"Q_wf (kW) = {Q_wf:.1f}")
 # print(f"Cp (J/kgC) = {Cp:.0f}")
  print("==========================================")

if len(sys.argv) < 2:
   print (f"Usage: python {sys.argv[0]} <sample time>")
else:
   ts = sys.argv[1]
   report ("htr", ts)
