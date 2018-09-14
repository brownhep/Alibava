import ROOT as R
from Methods_alibava import *
import time
import os.path

R.gROOT.ProcessLine('.L LandauFit.C+')

#filePrefix = sys.argv[1]    #ex: FTH200P_LTR10
#folder = sys.argv[2]	#Where data files are stored. ex: '/home/dtersegno/Raw_Data/LongTerm/FTH200P_R10/'
runMin = sys.argv[1]    #The minimum and maximum run numbers tolook at.
runMax = sys.argv[2]
runinfo = sys.argv[3]

#Get Run information from currentruninfo.txt
f = open(runinfo,"r")
filePrefix = f.readline()[:-1]
datafolder = f.readline()[:-1]
folder = f.readline()[:-1]
regline = f.readline()[:-1]
regions = regline.split(' ')
vline = f.readline()[:-1]
Volts = vline.split(' ')
print (filePrefix, folder, regions, Volts)
f.close()

#scale allows for conversion of ADC counts to another value. The 192k value is ~num electrons.
#scale = 192000/255
scale = 1

for runNum in range(int(runMin),int(runMax)+1):
	for volt in reversed(Volts):

		root_name = filePrefix + "_" + str(runNum)
		#Creates Root File for each voltage to store Plots 
		RFile_name = folder + root_name + '.root'
		if not os.path.isfile(RFile_name):
			root_name = filePrefix + "_T" + str(runNum)
			RFile_name = folder + root_name + '.root'
		print (RFile_name)

		for i in regions:				
			print ('Run #',runNum,', V: ', volt, ', Region: ', i)
			doSimpleLandau(RFile_name, int(i), start=10)  
			#doLandau(RFile_name, int(i), start=10)  



