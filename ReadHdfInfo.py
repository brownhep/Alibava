import numpy as np
import h5py
import sys

DataFile = sys.argv[1]
hdf=h5py.File(DataFile,'r')
ls=list(hdf.keys())
#print('list of datasets in the file:\n',ls)

data1=hdf.get('header')
data1_items = list(data1.items())

data2=hdf.get('scan')
data2_items = list(data2.items())

data3=hdf.get('events')
data3_items = list(data3.items())

print('the info of the header:\n',data1,'Groups in header:\n',data1_items,'\nthe info of the scan:\n',data2,'Groups in scan:\n',data2_items,'\nthe info of the events:\n',data3,'Groups in events:\n',data3_items,)

signal = data3.get('signal')
print('first signal data of the calibration run:\n',signal.value,signal.shape)

#start = data2.get('start')
#print('start data:\n', start.value,start.shape)

#end = data2.get('end')
#print('end data:\n',end.value, end.shape)

#value = data2.get('value')
#print('value data:\n', value.value, value.shape)

time = data3.get('time')
print('time data:\n', time.value[0],time.value, time.shape)

#clock = data3.get('clock')
#print('clock data:\n', clock.value, clock.shape) 



hdf.close()
