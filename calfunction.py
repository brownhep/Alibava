import numpy
import h5py
import sys
import ROOT as R
import time
import os.path
from datetime import datetime, date, time, timedelta
from array import array
from struct import *

def get_gains(file_name, Pedestal, RFile_name = '', T_Data=0, scale = 1):

   FilePrefix = file_name.replace('.bin','').replace('.dat','').replace('.hdf','')

   hcalp = R.TProfile2D("hcal","Calibration Histogram",288,0,288,32,0,32)
   hcalp.Sumw2()
   hcalp.SetOption("colz")
   hcaln = R.TProfile2D("hcal","Calibration Histogram",288,0,288,32,0,32)
   hcaln.Sumw2()
   hcaln.SetOption("colz")


   hTest = R.TH1F('Test','Test',100,-50,50)
   HGain = R.TH1D('Gain','Gain',256,0,256)
   HOffset = R.TH1D('Offset','Offset',256,0,256)
   HChi2 = R.TH1D('Chi2','Chi2',256,0,256)
   HBadCh = R.TH1D('BadCh','Bad Channels',256,0,256)
   HNsyCh = R.TH1D('NsyCh','Noisy Channels',256,0,256)
   HDisconCh = R.TH1D('DisconCh','Disconnected Channels',256,0,256)

   ical = array('i', [0])
   strip = array('i', [0])
   pcharge = array ('f', [0.])
   ncharge = array ('f', [0.])
   T_Chuck = array ('f', [0.])
   pgain = array ('f', [0.])
   offset = array ('f', [0.])
   chi2 = array ('f', [0.])

   RTree = R.TTree('rawcal', 'Raw Calibration Data')
   RTree.Branch('ical',ical,'ical/I')
   RTree.Branch('strip',strip,'strip/I')
   RTree.Branch('pcharge',pcharge,'pcharge/F')
   RTree.Branch('ncharge',ncharge,'ncharge/F')
   RTree.Branch('T_Chuck',T_Chuck,'T_Chuck/F')

   RTree2 = R.TTree('cal', 'Gain/Offset Calibration Data')
   RTree2.Branch('strip',strip,'strip/I')
   RTree2.Branch('pgain',pgain,'pgain/F')
   RTree2.Branch('offset',offset,'offset/F')
   RTree2.Branch('chi2',chi2,'chi2/F')

   Gains = []
   BadCh = []

   #f = open(file_name,"r")
   #Open Hdf5 file
   hdf = h5py.File(file_name,'r')
   
   #-----Reading Binary Header
   #Read Data header
   #while True:
   
   #     y1 = f.read(4)
   #     if not y1: break
   #     z1 = unpack("I",y1)
   #     if (i%100 == 0): print i, z1
   #     h_length = int(z1[0])
   #     y1 = f.read(4)
   #     z1 = unpack("I",y1)
        #if (i%100 == 0): print z1
   #y1 = f.read(18)
   #     z1 = unpack("3I3H",y1)
   #     #z2 = unpack("9H",y1)
   #     if (i%100 == 0): print "?2: ",z1, z2
        #print z1, z2
   #y1 = f.read(32)
   #y1 = f.read(256)
   #z1 = unpack("128H",y1)
   data3 = hdf.get('events')
   signal = data3.get('signal')
   event = 0
   for event_data in signal:
     icalstep = event/(signal.shape[0]/32)
     #event_data = signal.value[event]
     for j in range(256): 
        charge = event_data[j]-Pedestal.GetBinContent(j+1)
        if (j == 100) and (event%50==0): print (event,charge)
        if event%2 + j%2==1:
           hcalp.Fill(j,icalstep,charge)
        else:
           hcaln.Fill(j,icalstep,charge)
     hTest.Fill(event_data[j]-Pedestal.GetBinContent(j+1))
   #y1 = f.read(32)
   #y1 = f.read(256)
   #z1 = unpack("128H",y1)
   #for j in xrange(128): hcal.Fill(j+128,icalstep,z1[j])
     if (event%(signal.shape[0]/32) == 0): print('Data', event, signal.value[event])
     event += 1

   print(event, 'events')
   print(hcalp.GetBinContent(100,1))
   #Fill tree with averages from profile histogram
   for i in range(256):
        h1 = R.TH1D('FitProfile','Fit Profile',32,0,256)
        f1 = R.TF1()
        lf = R.TLinearFitter(1)
        for j in range(32):
                ical[0] = 8*j
                strip[0] = i
                pcharge[0] = hcalp.GetBinContent(i+1,j+1)
                ncharge[0] = hcaln.GetBinContent(i+1,j+1)
                print(hcalp.GetBinContent(i+1,j+1), Pedestal.GetBinContent(i+1))
                RTree.Fill()
                h1.Fill(ical[0],pcharge[0])
                #h2.Fill(ical[0],ncharge[0])
                #print "Charge", ical[0],charge[0]
                #lf.AddPoint(8.0*j,charge[0],1)
        h1.Fit("pol1")
        f1 = h1.GetFunction("pol1")
        strip[0] = i
        pgain[0] = f1.GetParameter(1)
        Gains.append(pgain[0])
        offset[0] = f1.GetParameter(0)
        chi2[0] = f1.GetChisquare()
        if chi2[0] > 12:
                HNsyCh.Fill(i,1)
                BadCh.append(i)
        if pgain[0] < 0.1:
                HBadCh.Fill(i,1)
                BadCh.append(i)
        if pgain[0] > 0.6 and chi2[0] > 4: 
                HDisconCh.Fill(i,1)
                BadCh.append(i)
        HGain.Fill(i,pgain[0])
        HOffset.Fill(i,offset[0])
        HChi2.Fill(i,chi2[0])
        RTree2.Fill()
        print (f1.GetParameter(0), f1.GetParameter(1), f1.GetChisquare())
        h1.Delete()


   RFile = R.TFile(FilePrefix + ".root","RECREATE")
   hcalp.Write()
   hcaln.Write()
   hTest.Write()
   HGain.Write()
   HOffset.Write()
   HChi2.Write()
   HNsyCh.Write()
   HBadCh.Write()
   HDisconCh.Write()
   RTree.Write()
   RTree2.Write()
   RFile.Close()

   return Gains, BadCh



