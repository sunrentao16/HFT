# read data -------------------------------------
import os 
import scipy.io
os.chdir(os.path.abspath('.')) # set working dictionary

#data = scipy.io.loadmat("Data/CSCO_20141103.mat")
data = scipy.io.loadmat("FB_20141103.mat")
#
dt = data["data"][0,0]
listName=['Event','SellVolume','SellPrice',
           'BuyVolume','BuyPrice']
#
for i in range(len(listName)):
    exec(listName[i] + ' = dt[' + str(i) + ']')

# add some comments
I changed something
