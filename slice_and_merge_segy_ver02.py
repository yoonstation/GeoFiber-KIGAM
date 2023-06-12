#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import h5py
from obspy import read, Trace, Stream, UTCDateTime
from obspy.core import AttribDict
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
from obspy.io.segy.core import _read_segy
import numpy as np
import sys
import matplotlib.pyplot as plt
from pkprocess import * 


# In[ ]:


import sys

# Check if command-line arguments are provided
if len(sys.argv) ==5:
    # The first argument (index 0) is the script name, so the actual arguments start from index 1
    #argument = sys.argv[1]
    t_start_first_file=sys.argv[1]
    t_end_last_file=sys.argv[2]
    tr_start = sys.argv[3]
    tr_end = sys.argv[4]
    #print("Argument:", argument)
    
    # Use the argument in your code
else:
    print("the number of argument is", len(sys.argv))
    print("Incorrect argument provided.")
    print("Please input the starttime of the file in the first file, and the endtime of the last file. + the fist trace and the last trace")
    print("example: python ~.py 5 5 2000 5000")
    quit()


# In[ ]:


# Set the very simple GUI environment 
from tkinter import filedialog as fd
input_files = fd.askopenfilenames()
input_files = sorted(input_files)
print(input_files)


# In[ ]:


ie=0
stream=Stream()
for input_file in input_files:
    print(input_file)
    if input_file[-4:] == "segy":
        print("Proceeding ...................................")
        tmp=Stream()
        tmp=_read_segy(input_file)
        
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
        dt_utc = UTCDateTime(int(yyyy), int(mm), int(dd), int(hh), int(mm2), int(ss), int(msec) )
        if ie == 0:
            slice_start = dt_utc + float(t_start_first_file)
            print(ie)
            print(slice_start)
        if ie == len(input_files)-1:
            slice_end = dt_utc + float(t_end_last_file)
            print(ie)
            print(slice_end)

        # Print the result 
        print(dt_utc)
        
        ii=0
        if ie == 0:
            stream = tmp

        # trace 단에서 처리해야함. 안 그러면 그냥 옆에 붙어버림.
        for trace in tmp:
            trace.stats.starttime = dt_utc
            if ie > 0:
                 stream[ii]=stream[ii]+trace
            ii=ii+1 
        ie=ie+1
        
# Time Slice 
stream.trim(slice_start, slice_end) # trim 으로 해야지 줄어듬. slice로 하면 잠깐 생겼다가 사라짐.

GT = stream[0].stats.starttime 
name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
output_file = name_flag + "."
print()

# Trace Slice 
tmp=Stream()
trace_start=int(tr_start)
trace_end=int(tr_end)

ii=0
for trace in stream:
    if ii > trace_start and ii < trace_end:
        tmp.append(trace)
    ii=ii+1
    
output_file=output_file+"segy"
tmp.write(output_file, format="SEGY", data_encoding=1, byteorder=sys.byteorder) 

