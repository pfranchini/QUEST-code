'''
Fridge5_AERO5_Run4_DM13_TWDAQ05_analysed.dat
============================================
['secs', 'date', 'time', 'TW26htnV', 'TW26width', 'TW26freq', 'TW26drive', 'TW26quadnV', 'TW26hwd', 'TW26Vxbgrd', 'TW26Vybgrd', 'TW27htnV', 'TW27width', 'TW27freq', 'TW27drive', 'TW27quadnV', 'TW27hwd', 'TW27Vxbgrd', 'TW27Vybgrd', 'TW28htnV', 'TW28width', 'TW28freq', 'TW28drive', 'TW28quadnV', 'TW28hwd', 'TW28Vxbgrd', 'TW28Vybgrd', 'TW29htnV', 'TW29width', 'TW29freq', 'TW29drive', 'TW29quadnV', 'TW29hwd', 'TW29Vxbgrd', 'TW29Vybgrd', 'TW30htnV', 'TW30width', 'TW30freq', 'TW30drive', 'TW30quadnV', 'TW30hwd', 'TW30Vxbgrd', 'TW30Vybgrd', 'index', 'datetime', 'TW26htV', 'TW27htV', 'TW28htV', 'TW29htV', 'TW30htV', 'TW26quadV', 'TW27quadV', 'TW28quadV', 'TW29quadV', 'TW30quadV', 'TW26VxbackgrV', 'TW27VxbackgrV', 'TW28VxbackgrV', 'TW29VxbackgrV', 'TW30VxbackgrV', 'TW26VybackgrV', 'TW27VybackgrV', 'TW28VybackgrV', 'TW29VybackgrV', 'TW30VybackgrV', 'TW26htVcorr', 'TW27htVcorr', 'TW28htVcorr', 'TW29htVcorr', 'TW30htVcorr', 'TW26driveV', 'TW27driveV', 'TW28driveV', 'TW29driveV', 'TW30driveV', 'TW26hwdSI', 'TW27hwdSI', 'TW28hwdSI', 'TW29hwdSI', 'TW30hwdSI', 'TW26velocity', 'TW27velocity', 'TW28velocity', 'TW29velocity', 'TW30velocity', 'TW26force', 'TW27force', 'TW28force', 'TW29force', 'TW30force', 'TW26power', 'TW27power', 'TW28power', 'TW29power', 'TW30power', 'TW26widththerm', 'TW27widththerm', 'TW28widththerm', 'TW29widththerm', 'TW30widththerm', 'TW26widththermAR', 'TW27widththermAR', 'TW28widththermAR', 'TW29widththermAR', 'TW30widththermAR', 'TW26temperature', 'TW27temperature', 'TW28temperature', 'TW29temperature', 'TW30temperature', 'TW26widthpar', 'TW27widthpar', 'TW28widthpar', 'TW29widthpar', 'TW30widthpar']
'''

import csv
import numpy as np
import matplotlib.pyplot as plt

with open('/home/franchini/Downloads/Fridge5_AERO5_Run4_DM13_TWDAQ05_analysed.dat', 'rt') as f:
  reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
  lineData = list()
  cols = next(reader)
  
  # Print names of variables
  print(cols)

  for col in cols:
    # Create a list in lineData for each column of data.
    lineData.append(list())

  for line in reader:
    for i in xrange(0, len(lineData)):
      # Copy the data from the line into the correct columns.
      lineData[i].append(line[i])

  data = dict()

  for i in xrange(0, len(cols)):
    # Create each key in the dict with the data in its column.
    data[cols[i]] = lineData[i]

print(data['TW27widthpar'])

#print(data['TW27temperature'])

plt.plot(data['secs'],data['TW27widthpar'])

#plt.plot(data['secs'],data['TW27temperature'])

plt.show()

