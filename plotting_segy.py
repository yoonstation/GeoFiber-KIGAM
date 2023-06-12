#!/usr/bin/env python
# coding: utf-8

# In[1]:


from obspy import read, Trace, Stream, UTCDateTime
from obspy.core import AttribDict
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
from obspy.io.segy.core import _read_segy
import numpy as np
import sys
import math
import matplotlib.pyplot as plt 
from pkprocess import * 


# In[2]:


# Set the very simple GUI environment 
from tkinter import filedialog as fd
input_files = fd.askopenfilenames()
input_files = sorted(input_files)
print(input_files)


# Plotting multiple files
ie=0
stream=Stream()
for input_file in input_files:
    print(input_file)
    if input_file[-4:] == "segy":
        print("Proceeding ...................................")
        tmp=Stream()
        tmp=_read_segy(input_file)
        
        ntraces, ns = int(tmp.count()), len(tmp[0])
        sampling_interval = tmp[0].stats.delta 
        print(ntraces, ns)
        print(tmp[0].stats.delta)
        SpatialSampling = 1 # dummy value because no spatial spaming information is inherited.
        xmin=tmp[0].stats.segy.trace_header.trace_sequence_number_within_line
        
        # trace 단에서 처리해야함. 안 그러면 그냥 옆에 붙어버림.
        ii=0
        data=np.zeros(shape=(ns,ntraces))
        for trace in tmp:
            data[:, ii]=trace.data
            ii=ii+1 
        
        # Get a time iniformation among the title strings 
        refchanum = -5
        time_string_size=-22
        date_time_string=input_file[time_string_size+refchanum:refchanum]
        print(date_time_string)
        
        # Replace the dot with a colon to match the ISO 8061 format
        # 구관이 명관임 그냥 뭔가 문제가 있는데 옛날 gcf manual ipynb 방식으로 처리
        print(date_time_string)
        yyyy=date_time_string[0:4]
        mm=date_time_string[4:6]
        dd=date_time_string[6:8]
        hh=date_time_string[9:11]
        mm2=date_time_string[11:13]
        ss=date_time_string[13:15]
        msec=date_time_string[16:22] #"00" # 6자리로 지정해줘야 문제가 없음... 백만스물하나 ... 
        print(yyyy, mm, dd, hh, mm2, ss, msec)
        
        # Create a UTCDateTime object 
        GT = UTCDateTime(int(yyyy), int(mm), int(dd), int(hh), int(mm2), int(ss), int(msec) )
        

        
       # Plot and Save Segys to pngs
        f2=0
    #    xmin=0
    #    xmax=(f2+(ntraces-1)*SpatialSampling)
        xmax=f2+(ntraces-1)
        tmax=(ns-1)*sampling_interval
        ratio="auto"
        clip = 90 
        
        plt.figure(figsize=[15,10])
        plt.xlabel("Trace number", fontsize="large")
        plt.ylabel("Time(s)", fontsize="large")
        
        data=perc_clip(data, clip)
        #plt.imshow(data, aspect=ratio, extent=[xmin, xmax, tmax, 0], cmap='bwr');
        plt.imshow(data, aspect=ratio, extent=[xmin, xmax, tmax, 0], cmap='bwr');
    
        
        clb=plt.colorbar()
        clb.ax.set_ylabel("Amplitude", loc='center')
        
        ie=ie+1
        name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
        output_file = name_flag + "."        
        output_file=output_file + "png"
        plt.savefig(output_file, dpi=600, bbox_inches='tight');    
        plt.show();


# In[ ]:




