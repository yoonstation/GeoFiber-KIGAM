#!/usr/bin/env python
# coding: utf-8

# In[1]:


import h5py
from obspy import read, Trace, Stream, UTCDateTime
from obspy.core import AttribDict
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
from obspy.io.segy.core import _read_segy
import numpy as np
import sys
import matplotlib.pyplot as plt
from pkprocess import * 


# In[2]:


# Set the very simple GUI environment 
from tkinter import filedialog as fd
input_files = fd.askopenfilenames()


# In[3]:


# Define the input and output file paths 
#input_file="0000653260_2023-02-06_01.24.58.96832.hdf5"
#input_file="./230222/0000653259_2023-02-06_01.24.57.46832.hdf5"
#output_file="example_hdf5.segy"


# In[4]:



def hdf5_to_segy (data, ntraces, sampling_interval, GT, output_file):
    # Create an ObsPy stream object 
    stream= Stream()
    
    timeflag= str("%03d" % GT.julday) + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)
    
    # Loop over each trace and add it to the stream
    for i in range(ntraces):
        trace = Trace(data=np.require(data[:,i], dtype=np.float32))
        trace.stats.starttime = GT
        trace.stats.delta = sampling_interval # segy는 0.06 ms가 최대 성긴 샘플링 & 3만 샘플 최대.
        # trace.stats.starttime = UTCDataTime ... file name stamp 이용해야 함. 원노트 obspy segy header 참조
        
        if not hasattr(trace.stats, 'segy.trace_header'):
            trace.stats.segy={}
        trace.stats.segy.trace_header  = SEGYTraceHeader()
        trace.stats.segy.trace_header.trace_sequence_number_within_line= i # SEQ number
        trace.stats.segy.trace_header.original_field_record_number= int(timeflag) # FFID
        #trace.stats.segy.trace_header.shotpoint_number= int(timeflag) # ShotpointNumber?
        trace.stats.segy.trace_header.energy_source_point_number= int(timeflag) # ShotpointNumber
        #trace.stats.segy.trace_header.ensemble_number= int(timeflag) # CMP 
        trace.stats.segy.trace_header.trace_number_within_the_original_field_record = i # Channel number 
        
        if i % 1000 == 0:
            print(i, "/", ntraces, 'appended')
        
        stream.append(trace)
    print("done ... ")
    
    stream.stats=AttribDict()
    stream.stats.textualfile_header = 'Textual Header'
    stream.stats.binary_file_header=SEGYBinaryFileHeader()
    stream.stats.binary_file_header.trace_sorting_code=5
    
    output_file=output_file + "segy"
    # write out to SEGY file
    stream.write(output_file, format="SEGY", data_encoding=1, byteorder=sys.byteorder) 
    
    return 

def outputs_pngs_bin(data, ntraces, SpatialSampling, ns, sampling_interval, output_file):
    # plot parameters 
    f2=0
    xmin=0
#    xmax=(f2+(ntraces-1)*SpatialSampling)
    xmax=f2+(ntraces-1)
    tmax=(ns-1)*sampling_interval
    ratio="auto"
    clip=90
    
    plt.figure(figsize=[15,10])
    plt.xlabel("Trace number", fontsize="large")
    plt.ylabel("Time(s)", fontsize="large")
    
    data = perc_clip(data, clip)
    plt.imshow(data, aspect=ratio, extent=[xmin, xmax, tmax, 0], cmap='bwr');
    
    clb=plt.colorbar()
    clb.ax.set_ylabel("Amplitude", loc='center')
    
    output_file=output_file + "png"
    plt.savefig(output_file, dpi=600, bbox_inches='tight');
    
    return 


# In[5]:


for input_file in input_files:
    print(input_file)
    if input_file[-4:] == "hdf5":
        print("Proceeding ...................................")
        # Read in HDF5 file and See the Metadata Structure
        with h5py.File(input_file, "r") as f:
            #data= f['NOISE1'][:] # TEST DATA for Display
            data= f['DAS'][:] 
            ns, ntraces = data.shape
            print("ns, ntraces", "=", ns, ntraces)
            for name in f:
                print("===========================================================================================")
                print(name)
                print("===========================================================================================")
                if name=="DAQ":
                    for item in f[name]:
                        print(item, "=", f[name][item][0])
                        if item == "RepetitionFrequency":
                            sampling_interval = 1./f[name][item][0]
                if name=="Interrogator":
                    for item in f[name]:
                        print(item, "=", f[name][item][0])
                if name=="Metadata":
                    for item in f[name]:
                        print(item, "=", f[name][item][0])
                        if item == 'Timestamp':
                            GT=UTCDateTime(f[name][item][0])
                if name=="ProcessingServer":
                    for item in f[name]:
                        print(item, "=", f[name][item][0])
                        if item == "SpatialSampling":
                            SpatialSampling = f[name][item][0]

        print("================================================================================================")
        #sampling_interval=2000*60*5/1200*sampling_interval
        #sampling_interval=sampling_interval*1/10 # For SEGY TEST ONLY
        print("sampling_interval = ", sampling_interval )
    else:
        print("Please select a proper file format")
        quit()
        
    # Writing SEGY in two ways ...[ns over 30,000 and else ...]
    limit_ns=30000
    if ns > limit_ns:
        nsplit=int(ns/limit_ns)
        
        isplit = 0
        while isplit < nsplit:
            if isplit == nsplit -1:
                tmpdata_end = data[ (isplit)*limit_ns :ns+1, :] # 파이썬은 [이상:미만] 임에 주의.
                print("isplit=", isplit, np.shape(tmpdata_end))
                ns_tmp, ntraces_tmp = np.shape(tmpdata_end)
                #output_file=input_file[:-4]+str("%02d" % isplit)+"."
                name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
                output_file = name_flag + "."
                print(output_file + "segy")
                hdf5_to_segy (tmpdata_end, ntraces_tmp, sampling_interval, GT, output_file)
                outputs_pngs_bin(tmpdata_end, ntraces_tmp, SpatialSampling, ns_tmp, sampling_interval, output_file)
                print(isplit+1, "/", nsplit, "is done totally")
                print("---------------------------------------------------------------------------")
                
            else:
                tmpdata = data[ (isplit)*limit_ns : (isplit+1)*limit_ns, :]
                print("isplot = ", isplit, np.shape(tmpdata))
                ns_tmp, ntraces_tmp = np.shape(tmpdata)
                #output_file=input_file[:-4]+str("%02d" % isplit)+"."
                name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
                output_file = name_flag + "."
                print(output_file + "segy")
                hdf5_to_segy (tmpdata, ntraces_tmp, sampling_interval, GT, output_file)
                outputs_pngs_bin(tmpdata, ntraces_tmp, SpatialSampling, ns_tmp, sampling_interval, output_file)
                print(isplit+1, "/", nsplit, "is done totally")
                print("---------------------------------------------------------------------------")

            isplit = isplit + 1 
            GT = GT + sampling_interval*(limit_ns)
            print("UTCDateTime  :", GT)
            
    else:
        #output_file=input_file[:-4]
        name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
        output_file = name_flag + "."
        print(output_file + "segy")
        hdf5_to_segy (data, ntraces, sampling_interval, GT, output_file)
        outputs_pngs_bin(data, ntraces, SpatialSampling, ns, sampling_interval, output_file)
        print("1/1", "is done totally")
        print("---------------------------------------------------------------------------")


# In[ ]:




