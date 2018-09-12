import ROOT as R
import sys
import time
import os.path
from datetime import datetime, date, time, timedelta
from array import array
from struct import *
import numpy
import h5py

#R.gROOT.ProcessLine('.L /user_data/agarabed/ARCS/src/C++/LandauFit.C+')
R.gROOT.ProcessLine('.L /Code/LandauFit.C+')

def get_tempdata(file_name, outfolder, filePrefix):

    time = array('f',[0.])
    dewp = array('f',[0.])
    Vpelt = array('f',[0.])
    Ipelt = array('f',[0.])
    chuckT = array('f',[0.])
    airT = array('f',[0.])
    RH = array('f',[0.])

    RTree = R.TTree('env','Environment Data')
    RTree.Branch('time',time,'time/F')
    RTree.Branch('dewp',dewp,'dewp/F')
    RTree.Branch('Vpelt',Vpelt,'Vpelt/F')
    RTree.Branch('Ipelt',Ipelt,'Ipelt/F')
    RTree.Branch('chuckT',chuckT,'chuckT/F')
    RTree.Branch('airT',airT,'airT/F')
    RTree.Branch('RH',RH,'RH/F')

    hchuck = R.TGraph()
    hdp = R.TGraph()
    hrh = R.TGraph()
    hair = R.TGraph()

    #Opens environmental data file
    print('Reading in environmental data from file: ', file_name)
    if os.path.isfile(file_name):
        env_file = open(file_name,'r')
        cnt = 0
        temp_data=[]
        t_data = [0 for x in xrange(7)]
        t0 = 0
 
        for line in env_file:
            env_data = line.split('\t')
            for i in xrange(7): t_data[i] = float(env_data[i])
            temp_data.append([])
            for i in xrange(7): temp_data[cnt].append(t_data[i])                
            if t0 == 0: t0 = float(env_data[0])/(60*60*24)
            time[0] = float(env_data[0])/(60*60*24) - t0
            dewp[0] = float(env_data[1])
            chuckT[0] = float(env_data[2])
            airT[0] = float(env_data[3])
            RH[0] = float(env_data[4])
            Vpelt[0] = float(env_data[5])
            Ipelt[0] = float(env_data[6])
            RTree.Fill()
            hchuck.SetPoint(cnt, time[0], chuckT[0])
            hdp.SetPoint(cnt, time[0], dewp[0])
            hrh.SetPoint(cnt, time[0], RH[0])
            hair.SetPoint(cnt, time[0], airT[0])
            cnt += 1

    #ginger: check the stuff following the out statement
    #I'm not sure if it should be part of the if statement or not
        if os.path.isfile(outfolder+filePrefix+"_Summary.root"): RFile = R.TFile(outfolder+filePrefix+"_Summary.root",'UPDATE')
        else:   
           RFile = R.TFile(outfolder+filePrefix+"_Summary.root",'CREATE')
        RTree.Write("",R.TObject.kOverwrite)
        hchuck.SetName("ChuckT")
        hchuck.SetTitle("Chuck Temp vs Time")
        hchuck.GetXaxis().SetTitle("Time (Days)")
        hchuck.GetYaxis().SetTitle("T Chuck (C)")
        hchuck.Write("",R.TObject.kOverwrite)
        hdp.SetName("DewPt")
        hdp.SetTitle("Dew Point vs Time")
        hdp.GetXaxis().SetTitle("Time (Days)")
        hdp.GetYaxis().SetTitle("Dew Point (C)")
        hdp.Write("",R.TObject.kOverwrite)
        hrh.SetName("RH")
        hrh.SetTitle("Humidity vs Time")
        hrh.GetXaxis().SetTitle("Time (Days)")
        hrh.GetYaxis().SetTitle("Humidity (%)")
        hrh.Write("",R.TObject.kOverwrite)
        hair.SetName("AirT")
        hair.SetTitle("Air Temp vs Time")
        hair.GetXaxis().SetTitle("Time (Days)")
        hair.GetYaxis().SetTitle("T Air (C)")
        hair.Write("",R.TObject.kOverwrite)
        RFile.Close()

    #ginger: I commented this part out. There doesn't seem to be an if statement that goes along with it
    else: 
        print('Environmental data file does not exist')
        temp_data = 0

    #print(temp_data)
    return temp_data


def get_temp(t_data, time, start=0):

    chuckT = [0,0]

    if t_data != 0:
        ti0 = t_data[0][0]
        Ti0 = t_data[0][2]
        for i in xrange(start,len(t_data)):
            #print('debug', i, t_data[i][0], time, t_data[i][2])
            if t_data[i][0] < time:
                ti0 = t_data[i][0]
                Ti0 = t_data[i][2]
            if t_data[i][0] < time:
                ti0 = t_data[i][0]
                Ti0 = t_data[i][2]
                #print('a',ti0,time,float(ti0)-float(time))
                continue
            else:
                ti = t_data[i][0]
                Ti = t_data[i][2]
                if ti == ti0: chuckT = [Ti0,0]
                else: chuckT = [Ti0 + (time-ti0)*(Ti-Ti0)/(ti-ti0),i-1]
                break

    return chuckT

# The funtion to read the data file to a root file ----------------------------------------
def get_preped(file_name, title, chip, RFile_name = '', T_Data=0, scale = 1, Pedestal2 = None):

    #Sets whether you want results in ADC or Electrons
    print('SCALE:',scale)
    
    #Opens data file with pedestal information
    hdf = h5py.File(file_name,'r')
    
# Now we have opened the file we want to read -----------------------------------------------

    #Creates the TProfile which stores the information about channel averages and noise
    Ped_Data = R.TProfile('Ped Data ' + title,'Ped Data ' + title,256,0,256)
    Ped_Data.Sumw2() #Calculates TRUE average
    Ped_Data.SetOption('colz')
    Ped_CM = R.TH1F('PrepedCM','Preped CM',256,0,256)
    Noise_Data = R.TProfile('Noise Data ' + title,'Noise Data ' + title,256,0,256)
    Noise_Data.Sumw2() #Calculates TRUE average
    Noise_Data.SetOption('colz')
    Ped_Data_Ali = R.TProfile('Ped Data Ali' + title,'Ped Data Ali' + title,256,0,256)
    Ped_Data_Ali.Sumw2() #Calculates TRUE average
    Ped_Data_Ali.SetOption('colz')
    Ped_CM = R.TH1F('PrepedCM','Preped CM',256,0,256)
    Noise_Data_Ali = R.TProfile('Noise Data Ali' + title,'Noise Data Ali' + title,256,0,256)
    Noise_Data_Ali.Sumw2() #Calculates TRUE average
    Noise_Data_Ali.SetOption('colz')


    Chip0_Data = R.TH2I('APV0','APV0',128,0,128,100,400,600)
    Chip1_Data = R.TH2I('APV1','APV1',128,0,128,100,400,600)
    Chip0_Data.SetOption('colz')
    Chip1_Data.SetOption('colz')
    Noise_Region = R.TProfile('RegNoise','Noise by Chip',2,0,2)

# Initializing everthing ----------------------------------------------------------------------

    #Creates TObject to store time
    PedStartTime = R.TVectorF(1)
    PedStartTemp = R.TVectorF(1)
    PedNsyEvts = R.TVectorF(1)

    stp = array('i',[0])
    ped = array('f',[0.])

    RTree = R.TTree('ppd','Ped_Data')
    RTree.Branch('stp',stp,'stp/I')
    RTree.Branch('ped',ped,'ped/F')

    #Reads in Data and stores in TProfile hi
    event_data = [] #stores data from 256 channels of 1 event
    event = 0 #Counts the number of events
    NoisyEvts = 0 #Counts the number of events

# Starting to read the data file ----------------------------------------------------------------

  # List of HDF5 data File
    ls = list(hdf.keys())
    print('List of datasets in the file:\n', ls)

  # Read Data header - Time
    
    #(y = f.read(16)
    #z = unpack("4I",y)
    #print (z)----------->This is for Binary
    data3=hdf.get('events')
    print('Info of the event:\n', data3, '\n')
    File_item = list(hdf.items())
    print('Items in the file:\n', File_item, '\n')    
 
    tmstamp =data3.get('clock')
    #runtype = int(z[2])
    #h_length = int(z[3])
    #print (h_length)
    print('Time data:\n', tmstamp.value, '\n Time data Length:', tmstamp.shape)

    #Read Pedestals
    #y = f.read(256*4)
    data1=hdf.get('header')
    print('Info of the header:\n', data1, '\n')
    pedAli = data1.get('pedestal')
    for i in range(256): Ped_Data_Ali.Fill(i,pedAli.value[0][i])
    print('Pedestal:\n',pedAli.value,'\n',pedAli.shape)

    #Read Noise
    noise = data1.get('noise')
    print('noise:\n',noise.value,'\n',noise.shape)
    for i in range(256): Noise_Data_Ali.Fill(i,noise.value[0][i])

    #Read Signal (full event)
    full_event = data3.get('signal')
    print('Signal:\n', len(full_event),full_event.value[0])
    
    #hdf.close()

# This is the part I haven't really understood-----------------------------
        #for j in xrange(128): hcal.Fill(j+128,icalstep,chip1[j])
    #full_event = chip0
    #full_event.extend(chip1)
    #if (i%100 == 0): print('Data', i, len(full_event))
    #i += 1
    event = 0
    for event_data in full_event: #full_event[:1000]:to run through 0-1000 events

##
#        event_data = full_event.value[event]
        if event % 1000 == 0: print (event,' ', event_data)

        event += 1
      


    #event += 1 #increments the event count
    #   if event % 4000 == 0: print(event/4) #Prints status

    #Converts to electrons if scale is not equal to 1
    #if scale != 1: event_data += [x*192000.0/full_event.value[-1] for x in full_event]   
    #else: event_data += full_event

        if chip==0: apv_event = event_data[0:127]
        else: apv_event = event_data[128:255]

#Fill TProfile with all 256 elements/channels of the event_data
#avgevt = sum(event_data)/float(len(event_data))
        avgevt = sum(apv_event)/float(len(apv_event))
        rmssum = 0
        for chan in range(len(apv_event)): rmssum += (apv_event[chan] - avgevt)*(apv_event[chan] - avgevt)                
        if (rmssum/float(len(apv_event))) > 200:
          NoisyEvts += 1
        else:        
          for chan in range(len(event_data)):
      #Fills the Pedestal with channel values
            if Pedestal2 == None: Ped_Data.Fill(chan, event_data[chan])
      #Option if you want to compare pedestals, not important
            else: Ped_Data.Fill(chan, event_data[chan]-Pedestal2.GetBinContent(chan+1))
            stp[0] = int(chan)
            ped[0] = event_data[chan]
            RTree.Fill()
    #  if (chan >= 0) and (chan<128): Chip0_Data.Fill(chan,event_data[chan])
    #  if (chan >= 128) and (chan<256): Chip1_Data.Fill(chan%128,event_data[chan])

    Ped_CM.Fill(rmssum/float(len(apv_event)))
    event_data = []

    #hdf.close()

    #Pedestal = Ped_Data.ProfileX(title)

    print('Noisy Events:',NoisyEvts)
    PedNsyEvts[0] = NoisyEvts

    for i in range(256):
        Noise_Region.Fill(i/128, Ped_Data.GetBinError(i+1))
        Noise_Data.Fill(i, Ped_Data.GetBinError(i+1))

    #Saves the Pedestal to the RootFile
    if RFile_name != '':
        RFile = R.TFile(RFile_name,'UPDATE')
        print(type(PedStartTime),type(tmstamp.value))
        PedStartTime[0] = tmstamp.value[0]
        PedStartTime.Write("PedStartTime")
        #PedStartTemp[0] = -99.0
        #PedStartTemp.Write("PedStartTemp")
        Ped_Data.Write()
        Ped_CM.Write()
        Noise_Data.Write()
        Ped_Data_Ali.Write()
        Noise_Data_Ali.Write()
        PedNsyEvts.Write("PedNoisyEvts")
        Chip0_Data.Write()
        Chip1_Data.Write()
        Noise_Region.Write()
        #RTree.Write()

        #Pedestal.Write()
        RFile.Close()

    #Returns Pedestal for future use
    hdf.close()
    return Ped_Data

def get_pedestal(file_name, region, PrePed, title, RFile_name, scale = 1):

    file = open(file_name, 'r')
    
    Ped_Data = R.TProfile('Ped Data ' + title,'Ped Data ' + title,256,0,256)
    Ped_Data.Sumw2()
    Ped_Data.SetOption('colz')
    
    Total_Noise = R.TH1F('Noise ' + title, 'Noise ' + title, 40*scale, -10*scale, 10*scale)
    Total_Noise.Sumw2()
    Total_Noise.SetOption('hist')

    stp = array('i',[0])
    ppd = array('f',[0.])
    ped = array('f',[0.])

    RTree = R.TTree('pd','Ped_Data')
    RTree.Branch('stp',stp,'stp/I')
    RTree.Branch('ped',ped,'ped/F')

    event_data = []
    event = 0
    for line in file:
        if line.count(',') > 3:
            event += 1
            if event % 4000 == 0: print(event/4)
            full_event = [float(x) for x in line.replace('\r\t','').split('\t')[-1].split(',')]
            if scale != 1: event_data += [x*192000.0/full_event[-1] for x in full_event[12:-1]]
            else: event_data += full_event[12:-1]
            if event % 4 == 0:
                subtracted_event = [event_data[i] - PrePed.GetBinContent(i+1) for i in range(256)]
                #print subtracted_event
                max_chan = subtracted_event.index(max(subtracted_event[32*(region+1):32*(region+2)]))
                for chan in range(len(event_data)):
                    if chan < max_chan - 1 or chan > max_chan + 1:
                        Ped_Data.Fill(chan, event_data[chan])
                        if chan >= 32 * (region +1) and chan < 32 * (region +2): Total_Noise.Fill(subtracted_event[chan])
                        stp[0] = chan
                        ped[0] = event_data[chan]
                        RTree.Fill()
                event_data = []
    file.close()

    if RFile_name != '':
        RFile = R.TFile(RFile_name,'UPDATE')
        Ped_Data.Write()
        Total_Noise.Write()
        RTree.Write()
        RFile.Close()

    return Ped_Data

def find_bad_chans(file_name, chip, Ped, title, RFile_name, scale = 1, extra = None):

    if extra == None: extra = 'Chip '+ str(chip)
    file = open(file_name, 'r')

    Sig_Data = R.TProfile(extra + ' Data',extra + ' Data',32, 0, 32)
    Sig_Data.Sumw2()
    Sig_Data.SetOption('colz')
    
    Sig_Data2 = R.TH2D(extra + ' Data2',extra + ' Data2',32, 0, 32, 150, -75.5*scale, 74.5*scale)
    Sig_Data2.Sumw2()
    Sig_Data2.SetOption('colz')
    
    good = R.TH1D('Good','Good',150, -75.5*scale,74.5*scale)
    bad = R.TH1D('Bad','Bad',150, -75.5*scale,74.5*scale)
    bad2 = R.TH1D('Bad2','Bad2',150, -75.5*scale,74.5*scale)
    
    good.Sumw2()
    bad.Sumw2()
    bad2.Sumw2()
    
    

    event_data = []
    event = 0
    for line in file:
        if line.count(',') > 3:
            event += 1
            if event % 4000 == 0: print(event/4)
            
            full_event = [float(x) for x in line.replace('\r\t','').split('\t')[-1].split(',')]
            if scale != 1: event_data += [x*192000.0/full_event[-1] for x in full_event[12:-1]]
            else: event_data += full_event[12:-1]
            
            if event % 4 == 0:
                subtracted_event = [event_data[i] - Ped.GetBinContent(i+1) for i in range(256)]
                reg_event = subtracted_event[128*chip:128*(chip+1)]
                for chan in range(len(reg_event)):
                    Sig_Data.Fill(chan, reg_event[chan])
                    Sig_Data2.Fill(chan, reg_event[chan])
                event_data = []

    dict = []
    bad_chans = []
    for i in range(32):
        if Sig_Data.GetBinContent(i+1) < (0)*scale:
            if i not in bad_chans: bad_chans.append(i)
            if i+1 not in bad_chans and i+1 < 32: bad_chans.append(i+1)
            if i-1 not in bad_chans and i-1 > -1: bad_chans.append(i-1)


    if RFile_name != '':
        RFile = R.TFile(RFile_name,'UPDATE')
        #Sig_Data.Write()
        #Sig_Data2.Write()
        #prof.Write()
        #Pedestal.Write()
        RFile.Close()
    print(bad_chans)
    return bad_chans

def get_signal(file_name, chip, Ped, title, RFile_name, bad_chans, RunNum, scale = 1, T_Array=0, stripstart=0, stripend=31, cut=0, extra = None):

    print ('strips', stripstart, stripend)
    if extra == None: extra = 'Chip'+str(chip)
    f = h5py.File(file_name, 'r')    

    runno = array('i', [0])
    eventno = array('i', [0])
    hit_strip = array('i', [0])
    clust_charge = array ('f', [0.])
    strip_charge = array ('f', [0.])
    tdctime = array ('f', [0.])
    #SR1_charge = array ('f', [0.])
    #SR2_charge = array ('f', [0.])
    #SL1_charge = array ('f', [0.])
    #SL2_charge = array ('f', [0.])
    R1 = array ('f', [0.])
    R2 = array ('f', [0.])
    L1 = array ('f', [0.])
    L2 = array ('f', [0.])
    eta = array ('f', [0.])
    CM_noise = array ('f', [0.])

    sig_noise = array ('f', [0.])
    soft_evt = array ('i', [0])
    T_Chuck = array ('f', [0.])


    RhdFile_name = RFile_name[:-5] + "_hd.root"
    print (RhdFile_name)
    RhdFile = R.TFile(RhdFile_name,'RECREATE')
    
    RTree = R.TTree('hits', 'Hit Data')
    RTree.Branch('runno',runno,'runno/I')
    RTree.Branch('eventno',eventno,'eventno/I')
    RTree.Branch('hit_strip',hit_strip,'hit_strip/I')
    RTree.Branch('clust_charge',clust_charge,'clust_charge/F')
    RTree.Branch('strip_charge',strip_charge,'strip_charge/F')
    RTree.Branch('time',tdctime,'tdctime/F')
    #RTree.Branch('SR1_charge',SR1_charge,'SR1_charge/F')
    #RTree.Branch('SR2_charge',SR2_charge,'SR2_charge/F')
    #RTree.Branch('SL1_charge',SL1_charge,'SL1_charge/F')
    #RTree.Branch('SL2_charge',SL2_charge,'SL2_charge/F')
    RTree.Branch('R1',R1,'R1/F')
    RTree.Branch('R2',R2,'R2/F')
    RTree.Branch('L1',L1,'L1/F')
    RTree.Branch('L2',L2,'L2/F')
    RTree.Branch('eta',eta,'eta/F')
    RTree.Branch('CM_noise',CM_noise,'CM_noise/F')
    RTree.Branch('sig_noise',sig_noise,'sig_noise/F')
    RTree.Branch('soft_evt',soft_evt,'soft_evt/I')
    RTree.Branch('T_Chuck',T_Chuck,'T_Chuck/F')

    if RFile_name != '':
        RFile = R.TFile(RFile_name,'UPDATE')

    CMs = R.TH1D('CMs '+extra,'CMs '+extra, 20, -5*scale, 5*scale)
    CMs.Sumw2()
    CMs.SetOption('hist')
    
    Sig_Data = R.TProfile(extra + ' Data',extra + ' Data',32, 0, 32)
    Sig_Data.Sumw2()
    Sig_Data.SetOption('colz')

    h_name_fine_bin = extra + ' Signal_' + '1000' + 'Bins'
    h_name_fine_bin_clust = extra + ' Signal_' + '1000' + 'Bins_clust'
    h_name_coarse_bin_clust = extra + ' Signal_' + '100' + 'Bins_clust'
    SignalHist1 = R.TH1D(h_name_fine_bin, h_name_fine_bin, 1000, 0, 100*scale)
    SignalHist1.Sumw2()
    SignalHist1.SetOption('hist')
    SignalHist2 = R.TH1D(h_name_fine_bin_clust, h_name_fine_bin_clust, 1000, 0, 100*scale)
    SignalHist2.Sumw2()
    SignalHist2.SetOption('hist')
    SignalHist3 = R.TH1D(h_name_coarse_bin_clust, h_name_coarse_bin_clust, 100, 0, 100*scale)
    SignalHist3.Sumw2()
    SignalHist3.SetOption('hist')
    EtaHist = R.TH1D("Eta","Eta", 120, -0.1, 1.1)
    EtaHist.Sumw2()
    EtaHist.SetOption('hist')

    #Creates TObject to store time
    StartTime = R.TVectorF(1)
    StartTemp = R.TVectorF(1)
    NoisyEvts = R.TVectorF(1)
    SoftEvts = R.TVectorF(1)
    NoiseHits = R.TVectorF(1)

    runno[0] = RunNum
    event_data = []
    StartTemp[0] = -99.0
    event = 0
    NbadChans = 0
    Nnoisehits = 0
    NsoftSignal = 0
    NnoisyEvts = 0
    NgoodEvts = 0
    Tindex = 0
    hit1 = hit2 = hit3 = 0

    print(bad_chans)

    #Read Data header
    data3=f.get('events')
    header=data3.get('header')
    print ('Header:\n', header.value,'\n', header.shape,'\n')
    tmstamp = data3.get('clock') 
    print ('Time:\n', tmstamp.value,'\n', tmstamp.shape,'\n')

    #Read Pedestals
    #y = f.read(256*4)
    data1=f.get('header')
    pedAli = data1.get('pedestal')
    print('Pedestals:\n',pedAli.value , '\n', pedAli.shape,'\n')

    #Read Noise
    noise = data1.get('noise')
    print('noise\n', noise.value, '\n', noise.shape,'\n')
   
    #Read Signal (full event)
    full_event = data3.get('signal')
    tdctime[0] = data3.get('time')
    print('Full Event:\n', full_event.value[0], '\n', full_event.shape,'\n')
 
    event = 0
    for event_data in full_event: #if run through 0-1000 events, type: full_event[:1000]
        event += 1
        print (event)
    tmstamp = 0
    tca = get_temp (T_Array, tmstamp, Tindex)
    T_Chuck[0] = tca[0]
    Tindex = tca[1]
    if StartTemp[0] < 90: StartTemp[0] = T_Chuck[0]

    if event % 4000 == 0: print(event/4) #print status
            
   # if scale != 1: event_data += [x*192000.0/full_event[-1] for x in full_event]
    #else: event_data += full_event
            
    eventno[0] = event

    subtracted_event = [event_data[i] - Ped.GetBinContent(i+1) for i in range(256)]
    apv_event = subtracted_event[chip*128:(chip+1)*128]
        #reg_event = subtracted_event[32*(region+1):32*(region+2)]
                
    for chan in range(len(apv_event)):
                Sig_Data.Fill(chan, apv_event[chan])

                #print max(subtracted_event), subtracted_event.index(max(subtracted_event))/32 - 1, max(reg_event), reg_event.index(max(reg_event))
                #raw_input('')
                
    eventSN = []
    for i in range(256):
                if Ped.GetBinError(i+1) == 0: eventSN.append(0)
                else: eventSN.append( subtracted_event[i]/float(Ped.GetBinError(i+1)*R.TMath.Sqrt(Ped.GetBinEntries(i+1))) )
            
            
    apvSN = eventSN[chip*128:(chip+1)*128]
        #regSN = eventSN[32*(region+1):32*(region+2)]

                #for i in xrange(len(regSN)):
                #        if regSN[i] > 2.75: hit3 +=1
                #        if regSN[i] > 2.5: hit2 += 1
                #        if regSN[i] > 2.25: hit1 += 1
                
    chan_nums = [i for i in range(128)]
                
    good_apvevt = [apv_event[i] for i in range(128) if i not in bad_chans]
    good_apvSN = [apvSN[i] for i in range(128) if i not in bad_chans]
    goodchan_nums = [chan_nums[i] for i in range(128) if i not in bad_chans]
                
                #Takes max signal to noise
                #good_max_reg_chan = good_regSN.index(max(good_regSN))
    good_max_apv_chan = good_apvevt.index(max(good_apvevt))
    max_apv_chan = goodchan_nums[good_max_apv_chan]
                
    if max_apv_chan in bad_chans: print(max_apv_chan)
    max_all_chan = chip*128 + max_apv_chan
    hit_strip[0] = max_apv_chan
    Noisehit = False
    if max_apv_chan < stripstart or max_apv_chan > stripend:
                Noisehit = True
                Nnoisehits += 1

    CM = sum(apv_event) - apv_event[max_apv_chan]
    count = 127.0

        #Loop over APV events and subtract off neighbors of the hit strip and bad/noisy channels
    for i in range(128):
        if abs(max_apv_chan - 1) == 1:
                        CM -= apv_event[i]
                        count -= 1
        if ((128*chip) + i) in bad_chans:
                        #print "removing bad_chan", i, apv_event[i]
                        CM -= apv_event[i]
                        count -= 1
                
        CM /= count
        CMs.Fill(CM)
        CM_noise[0] = CM
        NoisyEvt = False
        if abs(CM) > 5:        
                NoisyEvt = True
                NnoisyEvts+=1

        if max_apv_chan in bad_chans:
                event_data = []
                NbadChans+=1
                continue
            
                #if( apv_event[max_apv_chan] > 10): print apv_event[max_apv_chan], CM
                    
                #Fill unclustered histogram and tree brances
        seed_signal = apv_event[max_apv_chan] - CM
        total_signal = seed_signal
        strip_charge[0] = seed_signal

                #If hit is not on edges record entire 4 strip cluster
                #Possibly a bias if significant number of apv_evt[chan-1] = apv_evt[chan+1]
        wholecluster = False
        if (max_apv_chan > 0) and (max_apv_chan < 127): 
                if apv_event[max_apv_chan+1] > apv_event[max_apv_chan-1]:
                        if ((max_apv_chan + 2) < 128) and ((max_apv_chan - 1) > -1):
                                L1[0] = apv_event[max_apv_chan] - CM
                                R1[0] = apv_event[max_apv_chan+1] - CM
                                L2[0] = apv_event[max_apv_chan-1] - CM
                                R2[0] = apv_event[max_apv_chan+2] - CM
                                wholecluster = True
                else:
                        if ((max_apv_chan + 1) < 128) and ((max_apv_chan - 2) > -1):
                                L1[0] = apv_event[max_apv_chan-1] - CM
                                R1[0] = apv_event[max_apv_chan] - CM
                                L2[0] = apv_event[max_apv_chan-2] - CM
                                R2[0] = apv_event[max_apv_chan+1] - CM
                                wholecluster = True
                if apv_event[max_apv_chan+1] == apv_event[max_apv_chan-1]: print ('Left & Right neighbors equal charge!')


                #Do clustering and fill histograms. Change from 2 -> 1 times the noise to get counted in cluster
                #Now just add neighbors without doing a cut
        if (max_apv_chan + 1) < 128:
                total_signal += apv_event[max_apv_chan+1] - CM
                    #if (apv_event[max_apv_chan+1]-CM) > 1*Ped.GetBinError(max_all_chan + 1 + 1)*R.sqrt(Ped.GetBinEntries(max_all_chan + 1 + 1)): total_signal += apv_event[max_apv_chan+1] - CM
                    #SR1_charge[0] = apv_event[max_apv_chan+1] - CM
                    #if (max_reg_chan + 2) < 32:  SR2_charge[0] = apv_event[max_apv_chan+2] - CM
        if (max_apv_chan - 1) > -1:
                total_signal += apv_event[max_apv_chan-1] - CM
                    #if (apv_event[max_apv_chan-1]-CM) > 1*Ped.GetBinError(max_all_chan - 1 + 1)*R.sqrt(Ped.GetBinEntries(max_all_chan - 1 + 1)): total_signal += apv_event[max_apv_chan-1] - CM
                    #SL1_charge[0] = apv_event[max_apv_chan-1] - CM
                    #if (max_reg_chan - 2) > -1:  SL2_charge[0] = apv_event[max_apv_chan-2] - CM
                        
        clust_charge[0] = total_signal

                #Mark soft events (less than 3x noise) in tree (was previously 5x noise)
        soft_evt[0] = 0
                #if Ped.GetBinError(max_all_chan + 1)*R.sqrt(Ped.GetBinEntries(max_all_chan + 1)) == 0:
                #        print Ped.GetBinError(max_all_chan + 1), Ped.GetBinEntries(max_all_chan + 1), max_all_chan
                #        sig_noise[0] = 0
                #else: sig_noise[0] = float(reg_event[max_reg_chan])/(Ped.GetBinError(max_all_chan + 1)*R.sqrt(Ped.GetBinEntries(max_all_chan + 1)))
        sig_noise[0] = float(apv_event[max_apv_chan])/(Ped.GetBinError(max_all_chan + 1)*R.TMath.Sqrt(Ped.GetBinEntries(max_all_chan + 1)))
        if (apv_event[max_apv_chan]-CM) < cut*Ped.GetBinError(max_all_chan + 1)*R.TMath.Sqrt(Ped.GetBinEntries(max_all_chan + 1)):
                soft_evt[0] = 1
                #if (total_signal) < cut*Ped.GetBinError(max_all_chan + 1)*R.sqrt(Ped.GetBinEntries(max_all_chan + 1)):
                #    soft_evt[0] = 1

        eta[0] = -1
        if wholecluster and not(NoisyEvt) and not (Noisehit):
                R0 = R1[0]
                L0 = L1[0]
                #if R1[0] < 0:        R0 = 0
                #if L1[0] < 0:        L0 = 0
                if (R0 + L0) == 0:
                        print('zero alert',max_apv_chan,apv_event[max_apv_chan],Ped.GetBinError(max_all_chan + 1),Ped.GetBinEntries(max_all_chan + 1),CM,R1[0],L1[0])
                if soft_evt[0]==0: 
                        eta[0] = R0/(R0+L0)
                        EtaHist.Fill(eta[0])
                #print eta[0], R0, L0, R1, L1
                RTree.Fill()

        if (apv_event[max_apv_chan]-CM) < cut*Ped.GetBinError(max_all_chan + 1)*R.TMath.Sqrt(Ped.GetBinEntries(max_all_chan + 1)):
                NsoftSignal+=1
                event_data = []
                continue
                #if (total_signal) < cut*Ped.GetBinError(max_all_chan + 1)*R.sqrt(Ped.GetBinEntries(max_all_chan + 1)):
                #    NsoftSignal+=1
                #    event_data = []
                #    continue
                if NoisyEvt or Noisehit:
                    event_data = []
                    continue
                
        NgoodEvts+=1
        SignalHist1.Fill(seed_signal)
        SignalHist2.Fill(total_signal)
        SignalHist3.Fill(total_signal)

                
        event_data = []

    f.close()

    print('BAD EVENTS:',NbadChans)
    print('SOFT EVENTS:',NsoftSignal)
    print('NOISY EVENTS:',NnoisyEvts)
    print('NOISE HITS:',Nnoisehits)
    print('GOOD EVENTS:',NgoodEvts)
    print('Hits ', hit1, hit2, hit3)
    NoisyEvts[0] = NnoisyEvts
    SoftEvts[0] = NsoftSignal
    NoiseHits[0] = Nnoisehits

    if RFile_name != '':
        #RFile = R.TFile(RFile_name,'UPDATE')
        StartTime[0] = float(tmstamp)/(60*60*24) #Get time in days since Jan 1, 1904 + 40500 days
        StartTime.Write("StartTime",R.TObject.kOverwrite)
        StartTemp.Write("StartTemp",R.TObject.kOverwrite)
        NoisyEvts.Write("NoisyEvents",R.TObject.kOverwrite)
        SoftEvts.Write("SoftEvents",R.TObject.kOverwrite)
        NoiseHits.Write("NoiseHits",R.TObject.kOverwrite)
        CMs.Write("",R.TObject.kOverwrite)
        Sig_Data.Write("",R.TObject.kOverwrite)
        SignalHist1.Write("",R.TObject.kOverwrite)
        SignalHist2.Write("",R.TObject.kOverwrite)
        SignalHist3.Write("",R.TObject.kOverwrite)
        EtaHist.Write("",R.TObject.kOverwrite)
        #RTree = RTree.Write("",R.TObject.kOverwrite)
        RFile.Close()

    #RTree.Write("",R.TObject.kOverwrite)
    RhdFile.Write("",R.TObject.kOverwrite)
    RhdFile.Close()

def doSimpleLandau(RFile_name, chip, start=10, scale=1, extra=None): #start=85 was default in Methods, but in stream.py start=10 wass specified

        Landau_Chi2 = R.TVectorF(1)
        i = start


        if extra == None: extra = 'Chip'+str(chip)
        h_name_coarse_bin_clust = extra + ' Signal_' + '100' + 'Bins_clust'        
        hist_name = h_name_coarse_bin_clust

        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i)+')')
        RFile = R.TFile(RFile_name,'UPDATE')
        chi2 = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
        print(chi2, RFile.Get(hist_name).GetFunction(extra+'_Func').GetMaximumX(), i)
        Landau_Chi2[0] = chi2
        Landau_Chi2.Write("Chi2")
        RFile.Close()


def doLandau(RFile_name, chip, start=10, scale=1, extra=None): #start=85 was default in Methods, but in stream.py start=10 wass specified

        i = start

        if extra == None: extra = 'Chip'+str(chip)
        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i-1)+')')
        RFile = R.TFile(RFile_name,'READ')
        chi2_before = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
        RFile.Close()

        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i+1)+')')
        RFile = R.TFile(RFile_name,'READ')
        chi2_after = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
        RFile.Close()

        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i)+')')
        RFile = R.TFile(RFile_name,'READ')
        chi2 = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
        RFile.Close()
        print(chi2_before, chi2, chi2_after)

        if (chi2 > chi2_before and chi2 < chi2_after) or (chi2 > chi2_before and chi2 > chi2_after and chi2_before < chi2_after):
                while chi2 >= chi2_before or chi2 > chi2_after or chi2 > 500:
                        i -= 1.0
                        chi2_after = chi2
                        chi2 = chi2_before

                        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i)+')')
                        RFile = R.TFile(RFile_name,'READ')
                        chi2_before = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
                        print(chi2, RFile.Get(hist_name).GetFunction(extra+'_Func').GetMaximumX(), i)
                        RFile.Close()
                i += 1.0

        elif (chi2 < chi2_before and chi2 > chi2_after) or (chi2 > chi2_before and chi2 > chi2_after and chi2_before > chi2_after):
                while chi2 > chi2_before or chi2 >= chi2_after or chi2 > 500:
                        i += 1.0
                        chi2_before = chi2
                        chi2 = chi2_after
                        
                        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i)+')')
                        RFile = R.TFile(RFile_name,'READ')
                        chi2_after = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
                        print(chi2, RFile.Get(hist_name).GetFunction(extra+'_Func').GetMaximumX(), i)
                        RFile.Close()
                i -= 1.0

        R.gROOT.ProcessLine('LandauFit("' + RFile_name + '","' + hist_name + '","' + extra + '",'+str(scale)+','+str(i+1)+')')
        RFile = R.TFile(RFile_name,'READ')
        chi2_after = RFile.Get(hist_name).GetFunction(extra+'_Func').GetChisquare()
        print(chi2, RFile.Get(hist_name).GetFunction(extra+'_Func').GetMaximumX(), i)
        RFile.Close()
        

                
    
