import ROOT as R
from Methods_alibava import *
from calfunction import *
#from CalibrateSensors import *
import time
import numpy as np
import h5py 
#R.gROOT.ProcessLine('.L Code/LandauFit.C+')

#filePrefix = sys.argv[1]    #ex: FTH200P_LTR10
#folder = sys.argv[2]	#Where data files are stored. ex: '/home/dtersegno/Raw_Data/LongTerm/FTH200P_R10/'
runMin = sys.argv[1]    #The minimum and maximum run numbers tolook at.
runMax = sys.argv[2]

getEta = 0

#scale allows for conversion of ADC counts to another value. The 192k value is ~num electrons.
#scale = 192000/255
scale = 1

runno = array('i',[0])
rgn = array('i',[0])
strip = array('i',[0])
badch = array('i',[0])
pd = array('f',[0.])
noise = array('f',[0.])

RTree = R.TTree('strips','Strip Info')
RTree.Branch('run',runno,'runno/I')
RTree.Branch('rgn',rgn,'rgn/I')
RTree.Branch('strip',strip,'strip/I')
RTree.Branch('badch',badch,'badch/I')
RTree.Branch('pd',pd,'pd/F')
RTree.Branch('noise',noise,'noise/F')

	 
#Get Run information from currentruninfo.txt
f = open("currentruninfo.txt","r")
filePrefix = f.readline()[:-1]
folder = f.readline()[:-1]
outfolder = f.readline()[:-1]
regline = f.readline()[:-1]
chips = regline.split(' ')
vline = f.readline()[:-1]
Volts = vline.split(' ')
print(filePrefix, folder, chips, Volts)
f.close()

TFile_name = folder + filePrefix + '_Env.dat'
T_Array = get_tempdata(TFile_name, outfolder, filePrefix)

bad_chans = []

for runNum in range(int(runMin),int(runMax)+1):
	for volt in reversed(Volts):

		root_name = filePrefix + "_" + str(runNum)
		#Creates Root File for each voltage to store Plots 
		RFile_name = outfolder + root_name + '.root'
		RFile = R.TFile(RFile_name, 'RECREATE')
		RFile.Close()


		for i in chips:
			#Generates the Pre-Pedestal
			print('Prepedestal for ' + root_name +', ' + volt)
			#preped_name = folder + 'Raw_Data_'+root_name + '_Ped.dat'
			preped_name = folder + root_name + '_Ped.hdf'
			preped = get_preped(preped_name, 'PrePed', int(i), RFile_name, 1)
			print("Preped",preped)
			file_name = folder+root_name+'_Sr90.hdf'
			cal_name = folder+root_name+'_Cal.hdf'
			gain, bad_chans = get_gains(cal_name, preped) 
			if not os.path.isfile(file_name): file_name = folder + "Raw_Data_" + filePrefix + "_T" + str(runNum) + "_R" + i + ".hdf"
			#Generates Pedestal from Signal Data
			print ('Pedestal for ' + root_name +', ' + volt + ", R" + i)
			ped = get_pedestal(file_name, int(i), preped, 'R'+str(13-int(i)),RFile_name,scale)
			#Locates the bad Channels
			#bad_chans = []
			#bad_chans = find_bad_chans(file_name, int(i)i, ped, 'R'+str(13-int(i)),RFile_name,scale)
			#Calculates Signal + Fits Landau/Gaussian

			#Fill Tree of Strip Information
			for x in range(32):
				y = 32*(int(i)+1)+x
				runno[0] = runNum
				rgn[0] = int(i)
				strip[0] = y
				badch[0] = 0
				#if x in bad_chans:  badch[0] = 1
				pd[0] = ped.GetBinContent(y+1) #Was running off ped instead of preped
				noise[0] = ped.GetBinError(y+1)*R.TMath.Sqrt(preped.GetBinEntries(y+1))
				#print 'preped', x, preped.GetBinContent(y+1),preped.GetBinError(y+1)*R.sqrt(ped.GetBinEntries(y+1)), preped.GetBinEntries(y+1)
				#print 'ped   ', x, ped.GetBinContent(y+1),ped.GetBinError(y+1)*R.sqrt(ped.GetBinEntries(y+1)), ped.GetBinEntries(y+1)
				RTree.Fill()
				
			get_signal(file_name, int(i), ped, gain , 'R'+str(13-int(i)),RFile_name,bad_chans, runNum, scale, T_Array) # *****Try running off preped***** 


############ I comment this OUT!

	#		get_signal(file_name, int(i), preped, 'R'+str(13-int(i)),RFile_name,bad_chans, runNum, scale) # *****Try running off preped***** 





RStrip_Name = filePrefix + '_StripInfo.root'
RStripFile = R.TFile(RStrip_Name,'RECREATE')
RTree.Write()
RStripFile.Close()

