import sys
import ROOT
import numpy as np
import matplotlib as mpl
mpl.use('Agg')                    # Lets us plot without drawing the picture. Needed for running on cluster
import matplotlib.pyplot as plt
from matplotlib import cm


def counter(x):
  if ( (x % 100 == 0) and (x < 1000) ):
    print 'Entry', x
  elif ( (x % 1000 == 0) and (x < 10000) ):
    print 'Entry', x
  elif ( (x % 5000 == 0) and (x < 100000) ):
    print 'Entry', x
  elif ( x % 10000 == 0 ):
    print 'Entry', x

def loadrootfile(file):
  tfile = ROOT.TFile(file)
  tree = tfile.Get("deepid/DeePIDtree")
  tree._tfile = tfile # prevent tfile being deleted automatically
  return tree

file = sys.argv[1]

tree = loadrootfile(file)
entries = tree.GetEntriesFast()

print entries

wireLength = 500
timeLength = 500
timeMax = 4500
timeMin = 0
timeBinWidth = (timeMax-timeMin)/timeLength

for jentry in xrange(entries):
  if (jentry == 1000000000000000000):
    break
  counter(jentry)

  # get the next tree in the chain and verify
  ientry = tree.LoadTree(jentry)
  if ientry < 0:
    break

  # copy next entry into memory and verify
  nb = tree.GetEntry(jentry)
  if nb<=0:
    continue

  # use the values directly from the tree
  wires = tree.wireList
  times = tree.timeList
  tpcs = tree.TPCList
  adcs = tree.integralList
  planes = tree.planeList

  if ( len(wires) == 0 ):
    continue

  p0max = -1
  p0min = 1000000
  
  p1max = -1
  p1min = 1000000
  
  p2max = -1
  p2min = 1000000

  for i in range( len(wires[0]) ):
    wire = wires[0][i]
    plane = planes[0][i]
    if ( plane == 0 ):
      p0max = wire if wire > p0max else p0max
      p0min = wire if wire < p0min else p0min

    if ( plane == 1 ):
      p1max = wire if wire > p1max else p1max
      p1min = wire if wire < p1min else p1min
      
    if ( plane == 2 ):
      p2max = wire if wire > p2max else p2max
      p2min = wire if wire < p2min else p2min

  plane0 = np.zeros( ( timeLength , wireLength ) )
  plane1 = np.zeros( ( timeLength , wireLength ) )
  plane2 = np.zeros( ( timeLength , wireLength ) )

  for i in range( len(wires[0]) ):
    wire = wires[0][i]
    time = int( times[0][i]/timeBinWidth )
    plane = planes[0][i]
    adc = adcs[0][i]
    
    if ( plane == 0 ):
      wire -= p0min
      if ( wire >= 500 or time >= 500 ):
        continue
      plane0[time][wire] += adc

    if ( plane == 1 ):
      wire -= p1min
      if ( wire >= 500 or time >= 500 ):
        continue
      plane1[time][wire] += adc

    if ( plane == 2 ):
      wire -= p2min
      if ( wire >= 500 or time >= 500 ):
        continue
      plane2[time][wire] += adc

  wholeEvent = np.concatenate((plane0,plane1,plane2),axis=1)
#  rowsToRemove = 0
#  for index, row in enumerate(wholeEvent):
#    if ( max(row) == 0.0 ):
#      rowsToRemove += 1
#    elif( max(row) > 0.0 ):
#      break
#  wholeEvent = np.delete(wholeEvent, slice(0,rowsToRemove), axis=0)

  w = 15
  h = 5 #- rowsToRemove/100.0

  #Plotlocation = ("Plots/ele_" + str(jentry) + ".png") if (file == "DeePID_ele_merge.root") else ("Plots/muo_" + str(jentry) + ".png")
  Plotlocation = "Plots/" + sys.argv[2] + "_" + str(jentry) + ".png"
  fig = plt.figure(frameon=False)       # Remove the white space around the image. Ensures image only consists of data
  fig.set_size_inches(w,h)             # Set the size of the image, in inches
  ax = plt.Axes(fig, [0., 0., 1., 1.])  # Ensures image only consists of data
  ax.set_axis_off()                     # Ensures image only consists of data
  fig.add_axes(ax)                      # Add the object to the figure
  ax.imshow(wholeEvent,origin='lower',aspect='equal',cmap=cm.jet)  # Plot the numpy array, with the cm.jet colour map
  fig.savefig(Plotlocation,dpi=100)     # Save to the given location, at 100dpi. This gives us a 500x500 pixel image with the 5x5 inch size
  plt.close()








