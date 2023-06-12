#!/usr/bin/env python
# coding: utf-8

# In[1]:


from nptdms import TdmsFile
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


# In[3]:


# Define the input and output file paths 
#input_file="750m_gunpowder_50g_UTC+0900_DST0_20201210_095615.371.tdms"
#output_file="example_tdms.segy"


# In[4]:



def tdms_to_segy(data, ntraces,sampling_interval, GT, output_file):
    # Create an ObsPy stream object 
    stream= Stream()
    
    timeflag= str("%03d" % GT.julday) + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)
    
    # Loop over each trace and add it to the stream
    for i in range(ntraces):
        trace=Trace(data=np.require(data[:,i], dtype=np.float32))
        trace.stats.starttime = GT
        trace.stats.delta = sampling_interval # segy는 0.06 ms가 최대 성긴 샘플링 & 3만 샘플 최대.
        # trace.stats.starttime = UTCDataTime ... file name stamp 이용해야 함. 원노트 obspy segy header 참조
        
        if not hasattr(trace.stats, 'segy.trace_header'):
            trace.stats.segy={}
        trace.stats.segy.trace_header = SEGYTraceHeader()
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
    clip = 90 
    
    plt.figure(figsize=[15,10])
    plt.xlabel("Trace number", fontsize="large")
    plt.ylabel("Time(s)", fontsize="large")
    
    data=perc_clip(data, clip)
    plt.imshow(data, aspect=ratio, extent=[xmin, xmax, tmax, 0], cmap='bwr');
    
    clb=plt.colorbar()
    clb.ax.set_ylabel("Amplitude", loc='center')
    
    output_file=output_file + "png"
    plt.savefig(output_file, dpi=600, bbox_inches='tight');
    
    return 


# In[5]:


# Converting multiple files
for input_file in input_files:
    print(input_file)
    if input_file[-4:] == "tdms":
        print("Proceeding ...................................")
        # Read in TDMS file and See the Metadata Structure
        tdms_file=TdmsFile(input_file)
        metadata=TdmsFile.read_metadata(input_file)
        for item in metadata.properties:
            print(item, "=", metadata.properties[item])
            if item == "StreamTime[s]":
                if math.isinf(metadata.properties["StreamTime[s]"]) == True:
                    #print( "=========================================================================================")
                    ii=0
                    for group in tdms_file.groups():
                        group_name=group.name
                        #print(group_name)
                        ii=ii+1
                        iii=0
                        for channel in group.channels():
                            ns=int(len(channel[:]))
                            iii=iii+1
                            #print(ii, iii)
                            break 
                else:      
                    ns = int(metadata.properties["StreamTime[s]"] * metadata.properties["SamplingFrequency[Hz]"])
            if item == "SamplingFrequency[Hz]":
                sampling_interval = 1./metadata.properties["SamplingFrequency[Hz]"]
            if item == "SpatialResolution[m]":
                SpatialSampling = metadata.properties["SpatialResolution[m]"]
            if item ==  'GPSTimeStamp':
                GT=UTCDateTime(str(metadata.properties['GPSTimeStamp']))
        
        
        ntraces = int(metadata.properties["MeasureLength[m]"]/metadata.properties["SpatialResolution[m]"])
        timeflag= str("%03d" % GT.julday) + str("%02d" % GT.hour) + str(str("%02d" % GT.minute)) + str(str( "%02d" % GT.second))

        print(ns, ntraces, sampling_interval)
        
        data=np.zeros(shape=(ns, ntraces))
        
        for group in tdms_file.groups():
            group_name=group.name
            print(group_name)
            
            itrace=0
            for channel in group.channels():
                data[:,itrace]=channel[:]
                itrace=itrace+1     
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
                tmpdata_end = data[ (isplit)*limit_ns : ns+1, :] # 파이썬은 [이상:미만] 임에 주의.
                print("isplit=", isplit, np.shape(tmpdata_end))
                ns_tmp, ntraces_tmp = np.shape(tmpdata_end)
                #output_file=input_file[:-4]+str("%02d" % isplit)+"."
                name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
                output_file = name_flag + "."
                print(output_file + "segy")
                tdms_to_segy (tmpdata_end, ntraces_tmp, sampling_interval, GT, output_file)
                outputs_pngs_bin(tmpdata_end, ntraces_tmp, SpatialSampling, ns_tmp, sampling_interval, output_file)
                print(isplit, "/", nsplit, "is done totally")
                print("---------------------------------------------------------------------------")
                
            else:
                tmpdata = data[ (isplit)*limit_ns : (isplit+1)*limit_ns, :]
                print("isplot = ", isplit, np.shape(tmpdata))
                ns_tmp, ntraces_tmp = np.shape(tmpdata)
                #output_file=input_file[:-4]+str("%02d" % isplit)+"."
                name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
                output_file = name_flag + "."
                print(output_file + "segy")
                tdms_to_segy (tmpdata, ntraces_tmp, sampling_interval, GT, output_file)
                outputs_pngs_bin(tmpdata, ntraces_tmp, SpatialSampling, ns_tmp, sampling_interval, output_file)
                print(isplit, "/", nsplit, "is done totally")
                print("---------------------------------------------------------------------------")


            isplit = isplit + 1 
            GT = GT + sampling_interval*(limit_ns)
            print("UTCDateTime :", GT)
            
    else:
        #output_file=input_file[:-4]
        name_flag = str("%04d" % GT.year) + str("%02d" % GT.month) +  str("%02d" % GT.day) + "T" + str("%02d" % GT.hour) + str("%02d" % GT.minute) + str( "%02d" % GT.second)+"."+str("%06d" % GT.microsecond)
        output_file = name_flag + "."
        print(output_file + "segy")
        tdms_to_segy (data, ntraces, sampling_interval, GT, output_file)
        outputs_pngs_bin(data, ntraces, SpatialSampling, ns, sampling_interval, output_file)
        print("1/1", "is done totally")
        print("---------------------------------------------------------------------------")    


# In[ ]:






# In[ ]:





# In[ ]:




