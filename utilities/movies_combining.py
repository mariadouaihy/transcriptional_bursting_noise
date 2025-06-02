#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:55:45 2023

@author: mdouaihy
"""


import numpy as np
import pandas as pd
import os


def movies_combining(xlsfile,respath, filesList, fParam, extension):
    ###############################################

    fparam = np.load(fParam)

    FreqEchSimu = fparam['FreqEchSimu']
    FreqEchImg = fparam['FreqEchImg']
    TailleSeqMarq = fparam['TailleSeqMarq']
    TaillePostMarq = fparam['TaillePostMarq']
    Polym_speed = fparam['Polym_speed']
    EspaceInterPolyMin = fparam['EspaceInterPolyMin']
    FrameLen = fparam['FrameLen']
    DureeSignal = fparam['DureeSignal']
        
    
    tend = np.zeros((len(filesList)))
    tstart = np.zeros((len(filesList)))
    nfiles_l = np.ones((len(filesList)))
    true_FrameLen = np.zeros((len(filesList)))  # some movies may have a different time resolution
    
    
    ### this iteration is to fix the total number of nuclei and movie length
    ###d and also double checking the removal of 0 nuclei
    for iii in range(len(filesList)):
        fname = respath + 'result_' + filesList[iii].replace('_','') + '.npz'
        ffiles = np.load(fname)
        DataExp = ffiles['DataExp']
        DataPred = ffiles['DataPred']
        PosPred = ffiles['PosPred']
        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        
        n2 = DataExp.shape
    
        fname = xlsfile + filesList[iii] + extension

        ### loading time into mitosis
        rawData = pd.read_excel(fname)
        rawData.columns = [x.lower() for x in rawData.columns]
        time_into_mitosis = rawData['time'].dropna(how='all').to_numpy()
        
        
        if len(time_into_mitosis) ==0:
            true_FrameLen[iii] =  3.86 
            print('!!!!! framelen pre-assigned to 3.86')
            tstart[iii] = 0
            tend[iii] = n2[0]*3.86
        else:
            true_FrameLen[iii] = np.unique(np.round(np.diff(time_into_mitosis),3))
            tstart[iii] = time_into_mitosis[0]
            tend[iii] = time_into_mitosis[-1]
            
        nfiles_l[iii] = n2[1]
    
    
    ##############
    
    
    nnuclei = round(sum(nfiles_l)) # total number of nucleis
    cumsumnuclei = np.insert(np.cumsum(nfiles_l),0,0)
    cumsumnuclei = cumsumnuclei.astype(int)

    
    frame_num_combined = round(min(tend/true_FrameLen))
    
    
    DureeSimu = frame_num_combined*FrameLen ### film duration in s
    DureeAnalysee = DureeSimu + DureeSignal;
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed))

    DataExp_combined = np.zeros((frame_num_combined, nnuclei))
    DataPred_combined= np.zeros((frame_num_combined, nnuclei))
    PosPred_combined= np.zeros((num_possible_poly, nnuclei))
    nn_combined = PosPred_combined.shape
    tmax_combined = np.ones((nnuclei,))
 
    
    #### combining movies according to their time into mitosis
    for iii in range(len(filesList)):

        fname = respath + 'result_' + filesList[iii].replace('_','') + '.npz'
        ffiles = np.load(fname)
        DataExp = ffiles['DataExp']
        DataPred = ffiles['DataPred']
        PosPred = ffiles['PosPred']
        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        DataPred = np.delete(DataPred,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        DataPred = np.delete(DataPred, check_nan_columns, 1)
        PosPred = np.delete(PosPred, check_nan_columns, 1)
        
        n2 = DataExp.shape
        n3 = PosPred.shape   
        ############## concatenating nuclei according to time into mitosis
        if tstart[iii] == 0:

            DataExp_combined[:n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:frame_num_combined,:]
            DataPred_combined[:n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataPred[:frame_num_combined,:]
            t0_posPred = 0
            PosPred_combined[t0_posPred: t0_posPred+n3[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = PosPred[:num_possible_poly,:]

            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]

        else:

            DataExp_combined[round(tstart[iii]/FrameLen): round(tstart[iii]/FrameLen)+n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:frame_num_combined-round(tstart[iii]/FrameLen), :]
            DataPred_combined[round(tstart[iii]/FrameLen): round(tstart[iii]/FrameLen)+n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataPred[:frame_num_combined-round(tstart[iii]/FrameLen), :]
            
            t0_posPred = round(tstart[iii]*FreqEchSimu)
            
            if (t0_posPred+n3[0]-num_possible_poly) <=0:
                PosPred_combined[t0_posPred: t0_posPred+n3[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = PosPred[:num_possible_poly-t0_posPred, :]
            else:
                PosPred_combined[t0_posPred: t0_posPred+n3[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = PosPred[:-(t0_posPred+n3[0]-num_possible_poly),:]
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]
        
    
        if round(tstart[iii]/FrameLen)+n2[0] != frame_num_combined:
            DataExp_combined[round(tstart[iii]/FrameLen)+n2[0]:,cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan
            DataPred_combined[round(tstart[iii]/FrameLen)+n2[0]:, cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan
            PosPred_combined[t0_posPred+n3[0]-1:, cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan


    T0_combined = np.zeros((nn_combined[1],))  ##### will contain the start of the analyzed region

    PosPred_combined_corrected = PosPred_combined.copy()
    ##### correct PosPred: eliminate positions not contributing to intensity
    ##### compute T0_combined
    for i in range(nn_combined[1]): ### for all cells
        pospol =  np.where(PosPred_combined[:, i] == 1)[0]
        times = pospol / FreqEchSimu   ### starting times of polymerases in seconds
        ################################## are at times -  (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed)
        ################################## so seq start at times -  (TailleSeqMarq+TaillePostMarq)/Polym_speed)
        ind = np.where(times -  (TailleSeqMarq+TaillePostMarq)/Polym_speed > tmax_combined[i])[0]  # ### positions that have no influence
        PosPred_combined_corrected[pospol[ind], i] = 0
        max_intensity = max(DataPred_combined[:, i])
        ihit = np.where(DataPred_combined[:, i] > max_intensity / 5)[0]
        if len(ihit) != 0:
            ihi_min =  min(ihit)
            T0_combined[i] = (ihi_min-1)/FreqEchImg  #### T0_combined
        else:
            T0_combined[i] = tmax_combined[i]

    return [DataExp_combined, DataPred_combined, PosPred_combined, PosPred_combined_corrected, T0_combined, tmax_combined, FrameLen]










def movies_combining_ponxkini(xlspath, resultpath, fParam, extension, tfinal =1000, Tstart = 0 ):
 
    content                                                      = np.load(fParam)
    FreqEchSimu                                                  = content['FreqEchSimu']
    FreqEchImg                                                   = content['FreqEchImg']
    DureeSignal                                                  = content['DureeSignal']
    
    """  combine  movies """
      
    files                                                        = np.array(os.listdir(xlspath))
    files                                                        = list(map(lambda x: x.replace(extension,'') ,files))
    
    [DataExp, DataPred, PosPred_not_corrected, PosPred, T0_movie, tmax_movie, FrameLen] = movies_combining(xlspath, resultpath, files, fParam, extension)

    DataExp                                                      =  DataExp[int(Tstart*60/FrameLen):min(int(tfinal*60/FrameLen), DataExp.shape[0]),:]
    DataPred                                                     = DataPred[int(Tstart*60/FrameLen):min(int(tfinal*60/FrameLen), DataPred.shape[0]),:]
    PosPred                                                      = PosPred[int(Tstart*60*FreqEchSimu):min(int((DureeSignal + tfinal*60)*FreqEchSimu), PosPred.shape[0]),:]
    
    check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
    DataExp = np.delete(DataExp, check_nan_columns, 1)
    DataPred = np.delete(DataPred, check_nan_columns, 1)
    PosPred = np.delete(PosPred, check_nan_columns, 1)
    
    
    T0                                                           = np.zeros(DataExp.shape[1])
    tmax                                                         = np.zeros(DataExp.shape[1])
    for n_i in range(DataExp.shape[1]):
        ihit        = np.where(DataPred[:,n_i]> np.nanmax(DataPred[:,n_i])/5)[0]
        if len(ihit)== 0:
            if np.sum(DataExp[:,n_i]) ==0:
                print('here')
                T0[n_i]                   = np.nan
                tmax[n_i]                 = np.nan
            else:
                T0[n_i]                   = np.where(DataExp[:,n_i]!=0)[0][0] #min(tfinal*60, DataExp.shape[0]*FrameLen)
                tmax[n_i]                 = np.where(DataExp[:,n_i]!=0)[0][-1]
        else:
            T0[n_i]       = (min(ihit)+1)/FreqEchImg
            tmax[n_i]     = (max(ihit)+1)/FreqEchImg
    

    return T0, tmax, DataExp, DataPred, PosPred

    

def data_reshape(xlsFile, resultPath_i, fParam, tfinal = 10000, Tstart=0):
    """ load parameters """
    content           = np.load(fParam)
    DureeSignal       = content['DureeSignal']
    Polym_speed       = content['Polym_speed']
    EspaceInterPolyMin= content['EspaceInterPolyMin']
    FreqEchSimu       = content['FreqEchSimu']
    FreqEchImg        = content['FreqEchImg']
    TailleSeqMarq     = content['TailleSeqMarq']
    TaillePostMarq    = content['TaillePostMarq']
    
    """ load data """
    content           = np.load(resultPath_i)
    DataExp           = content['DataExp']
    DataPred          = content['DataPred']
    PosPred           = content['PosPred']



    """ loading time into mitosis""" 
    rawData           = pd.read_excel(xlsFile)
    rawData.columns   = [x.lower() for x in rawData.columns]
    time_into_mitosis = rawData['time'].dropna(how='all').to_numpy()
        
    if len(time_into_mitosis) == 0:
        print('!!!!! framelen pre-assigned to 3.86')
        tstart                  = 0
        tend                    = DataExp.shape[0]*3.86
        FrameLen                = 3.86
    else:
        tstart                  = time_into_mitosis[0]
        tend                    = time_into_mitosis[-1]
        FrameLen                = np.unique(np.round(np.diff(time_into_mitosis),3))[0]
    
    """ reshape the data to start from mitosis """
    frame_num         = round(tend/FrameLen)-1    
    DureeSimu         = frame_num*FrameLen ### film duration in s
    DureeAnalysee     = DureeSimu + DureeSignal
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed))
      
    n3                = PosPred.shape

    if tstart!=0:
        dataExp           = np.zeros((frame_num+1, DataExp.shape[1]))
        dataPred          = np.zeros((frame_num+1, DataExp.shape[1]))
        posPred           = np.zeros((num_possible_poly, DataExp.shape[1]))
    
        dataExp[int(np.round(tstart/FrameLen)): round(tstart/FrameLen)+frame_num, :]  = DataExp[:min((frame_num-round(tstart/FrameLen))+1, DataExp.shape[0]),:]
        dataPred[int(np.round(tstart/FrameLen)): round(tstart/FrameLen)+frame_num, :] = DataPred[:min((frame_num-round(tstart/FrameLen))+1, DataPred.shape[0]),:]
        
        t0_posPred                                                                    = round(tstart*FreqEchSimu)
        if (t0_posPred+n3[0]-num_possible_poly) <=0:
            posPred[t0_posPred: t0_posPred+n3[0], :] = PosPred[:num_possible_poly-t0_posPred, :]
        else:
            posPred[t0_posPred: t0_posPred+n3[0], :] = PosPred[:-(t0_posPred+n3[0]-num_possible_poly),:]

    
        DataExp           = dataExp.copy()
        DataPred          = dataPred.copy()
        PosPred           = posPred.copy()

    """ stop the data at either the minute pre-given or till the end of the movie """
    DataExp           = DataExp[int(Tstart*60/FrameLen):min(int(tfinal*60/FrameLen), DataExp.shape[0]),:]
    DataPred          = DataPred[int(Tstart*60/FrameLen):min(int(tfinal*60/FrameLen), DataPred.shape[0]),:]
    PosPred           = PosPred[int(Tstart*60*FreqEchSimu):min(int((DureeSignal + tfinal*60)*FreqEchSimu), PosPred.shape[0]),:]
    
    """ remove 0 nuclei """
    indx              = np.where(np.sum(DataPred, axis=0)!=0)[0]    
    DataExp           = DataExp[:,indx]
    DataPred          = DataPred[:,indx]
    PosPred           = PosPred[:,indx]
    
    """ get the time point where the signal reached 20 % of maximum intensity """
    T0                = np.zeros(DataExp.shape[1])
    tmax              = np.zeros(DataExp.shape[1])
    for i in range(DataExp.shape[1]): ### for all cells
        pospol                 =  np.where(PosPred[:,i] == 1)[0] 
        times                  = pospol / FreqEchSimu  
        ind                    = np.where( times -  (TailleSeqMarq+TaillePostMarq)/Polym_speed > tend )[0] #### positions that have no influence
        PosPred[pospol[ind],i] = 0 
        ihit                   = np.where(DataPred[:,i]> np.nanmax(DataPred[:,i])/5)[0]
        if len(ihit)!=0:
#            ihi_min                 =  min(ihit)
            T0[i]                   = (min(ihit)+1)/FreqEchImg;
            tmax[i]                 = (max(ihit)+1)/FreqEchImg;
        else:
            print(np.where(DataPred[:,i]!=0))
            T0[i]                   = np.where(DataPred[:,i]!=0)[0] #min(tfinal*60, DataExp.shape[0]*FrameLen)
            tmax[i]                 = np.where(DataPred[:,i]!=0)[-1]

    return T0, tmax, DataExp, DataPred, PosPred





def data_reshape_raw(xlsFile, tfinal = 1000, Tstart=0):
    """ load data """
    # First read the column names to determine available columns
    df_columns = pd.read_excel(xlsFile, nrows=0, engine='openpyxl')
    n_columns = len(df_columns.columns)
    # Use only columns from index 4 up to the actual number of columns
    usecols = np.arange(4, min(n_columns, 300))
    
    DataExp = pd.read_excel(xlsFile, usecols=usecols, skiprows=range(0,1), engine='openpyxl').to_numpy()
    
    ### removing zeros nuclei
    idx_zeros = np.where(np.sum(DataExp,axis=0)==0)[0] #avoid zero cells
    DataExp = np.delete(DataExp, idx_zeros, axis=1)

    
    check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
    DataExp = np.delete(DataExp,check_nan_rows,0)
    
    check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
    DataExp = np.delete(DataExp, check_nan_columns, 1)
        


    """ loading time into mitosis""" 
    rawData           = pd.read_excel(xlsFile)
    rawData.columns   = [x.lower() for x in rawData.columns]
    time_into_mitosis = rawData['time'].dropna(how='all').to_numpy()
        
    if len(time_into_mitosis) == 0:
        print('!!!!! framelen pre-assigned to 3.86')
        tstart                  = 0
        tend                    = DataExp.shape[0]*3.86
        FrameLen                = 3.86
    else:
        tstart                  = time_into_mitosis[0]
        tend                    = time_into_mitosis[-1]
        FrameLen                = np.unique(np.round(np.diff(time_into_mitosis),3))[0]
    
    """ reshape the data to start from mitosis """
    frame_num         = round(tend/FrameLen)-1    
    
    if tstart!=0:
        dataExp           = np.zeros((frame_num+1, DataExp.shape[1]))
        dataExp[int(np.round(tstart/FrameLen)): round(tstart/FrameLen)+frame_num, :]  = DataExp[:min((frame_num-round(tstart/FrameLen))+1, DataExp.shape[0]),:]
        
        DataExp           = dataExp.copy()
    
    """ stop the data at either the minute pre-given or till the end of the movie """
    DataExp           = DataExp[int(Tstart*60/FrameLen):min(int(tfinal*60/FrameLen), DataExp.shape[0]),:]
    
    
    """ get the time point where the signal reached 20 % of maximum intensity """
    T0                = np.zeros(DataExp.shape[1])
    tmax              = np.zeros(DataExp.shape[1])
    for i in range(DataExp.shape[1]): ### for all cells 
        ihit                   = np.where(DataExp[:,i]> np.nanmax(DataExp[:,i])/5)[0]
        if len(ihit)!=0:
            T0[i]                   = (min(ihit)+1)*FrameLen
            tmax[i]                 = (max(ihit)+1)*FrameLen
        else:
            print(DataExp[:,i])
            T0[i]                   = np.where(DataExp[:,i]!=0)[0] #min(tfinal*60, DataExp.shape[0]*FrameLen)
            tmax[i]                 = np.where(DataExp[:,i]!=0)[-1]

    return T0, tmax, DataExp


def movies_combining_rawData(xlsPath, calibrtion, extension):


    """ extracting files names """
    file_name_list = np.array(os.listdir(xlsPath))     # list of the data
    filesList = list(map(lambda x: x.replace(extension, ''), file_name_list))


    ###############################################

    tend = np.zeros((len(file_name_list)))
    tstart = np.zeros((len(file_name_list)))
    true_FrameLen = np.zeros((len(file_name_list)))   # some movies may have a different time resolution
    nfiles_l = np.ones((len(file_name_list)))

    ### this iteration is to fix the total number of nuclei and movie length
    ###d and also double checking the removal of 0 nuclei 
    for iii in range(len(file_name_list)):
        fname = xlsPath +  filesList[iii] + extension

        # First read the column names to determine available columns
        df_columns = pd.read_excel(fname, nrows=0, engine='openpyxl')
        n_columns = len(df_columns.columns)
        # Use only columns from index 4 up to the actual number of columns
        usecols = np.arange(4, min(n_columns, 300))
        
        DataExp = pd.read_excel(fname, usecols=usecols, skiprows=range(0,1), engine='openpyxl').to_numpy()
        
        ### removing zeros nuclei
        idx_zeros = np.where(np.sum(DataExp,axis=0)==0)[0] #avoid zero cells
        DataExp = np.delete(DataExp, idx_zeros, axis=1)

        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1) == DataExp.shape[1])[0]
        DataExp = np.delete(DataExp, check_nan_rows, 0)

        check_nan_columns = np.where(np.sum(DataExp == 0, axis=0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)

        DataExp = DataExp*calibrtion

        n2 = DataExp.shape

        # ## loading time into mitosis
        rawData = pd.read_excel(fname, engine='openpyxl')
        rawData.columns = [x.lower() for x in rawData.columns]
        time_into_mitosis = rawData['time'].dropna(how='all').to_numpy()

        if len(time_into_mitosis) == 0:
            true_FrameLen[iii] =  3.86
            print('!!!!! framelen pre-assigned to 3.86')
            tstart[iii] = 0
            tend[iii] = n2[0]*3.86
        else:
            true_FrameLen[iii] = np.unique(np.round(np.diff(time_into_mitosis),3))
            tstart[iii] = time_into_mitosis[0]
            tend[iii] = time_into_mitosis[-1]
        
        nfiles_l[iii] = n2[1]
    
    
    ##############
    
    
    nnuclei = round(sum(nfiles_l)) # total number of nucleis
    cumsumnuclei = np.insert(np.cumsum(nfiles_l),0,0)
    cumsumnuclei = cumsumnuclei.astype(int)
    
    
    frame_num_combined = round(min(tend/true_FrameLen))-1
    
    dataExp = np.zeros((frame_num_combined+1, nnuclei))
    tmax_combined = np.ones((nnuclei,))
    
    FrameLen = np.unique(true_FrameLen)[0]
    
    
    #### combining movies according to their time into mitosis
    for iii in range(len(filesList)):
        fname = xlsPath +  filesList[iii] + extension
        
        # First read the column names to determine available columns
        df_columns = pd.read_excel(fname, nrows=0, engine='openpyxl')
        n_columns = len(df_columns.columns)
        # Use only columns from index 4 up to the actual number of columns
        usecols = np.arange(4, min(n_columns, 300))
        
        DataExp = pd.read_excel(fname, usecols=usecols, skiprows=range(0,1), engine='openpyxl').to_numpy()
        
        ### removing zeros nuclei
        idx_zeros = np.where(np.sum(DataExp,axis=0)==0)[0] #avoid zero cells
        DataExp = np.delete(DataExp, idx_zeros, axis=1)

        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        
        n2 = DataExp.shape 
        
    ############## concatenating nuclei according to time into mitosis
        if tstart[iii]==0:
            
            dataExp[:frame_num_combined, cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:min(frame_num_combined, DataExp.shape[0]),:]
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]
        else:
#            shape
            dataExp[int(np.round(tstart[iii]/FrameLen)): round(tstart[iii]/FrameLen)+frame_num_combined, 
                    cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:min((frame_num_combined-round(tstart[iii]/FrameLen))+1, DataExp.shape[0]),:]
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]
    
        if np.round(tstart[iii]/FrameLen)+n2[0] != frame_num_combined:
            dataExp[int(np.round(tstart[iii]/FrameLen))+n2[0]:,cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan

            
    return [FrameLen, dataExp, tmax_combined, tstart]






def allign_data_Activation(data):
    #################### combine data according to time to reach 20% of full intensity
    data_aligned    = []
    for nuclei_i in range(data.shape[1]):
        ihit                = np.where(data[:,nuclei_i] > np.nanmax(data[:,nuclei_i])/5 )[0]
        if len(ihit)!=0:
            data_aligned.append(data[ihit[0]:,nuclei_i])
        else:
            data_aligned.append(data[:,nuclei_i])
            
    max_length      = max(len(row) for row in data_aligned)
    data_aligned    = np.array([list(row) + [np.nan] * (max_length - len(row)) for row in data_aligned], dtype=float).T        
    return data_aligned
            
def data_average_in_TW(data, TW_bp):
    sz              = data.shape
    TW_bp           = int(TW_bp)
    if len(sz) == 1:
        print('here')
        data_in_TW          = data[:sz[0] // TW_bp * TW_bp].reshape(sz[0] // TW_bp, TW_bp)
    else:        
        data_in_TW          = data[:sz[0] // TW_bp * TW_bp].reshape(sz[0] // TW_bp, TW_bp, sz[1])
    data_in_TW      = np.nanmean(data_in_TW, axis=1)            
    return data_in_TW


def clean_movie_name(movie_name_i):
    movie_name_i = (movie_name_i.replace('/', ' ')).replace('_',' ')
    movie_name_i = movie_name_i.replace('CalibratedTraces', '')
    movie_name_i = movie_name_i.replace('Calibrated', '')
    movie_name_i = movie_name_i.replace('nuclei', '')
    movie_name_i = movie_name_i.replace('time', '')
    movie_name_i = movie_name_i.replace('sorted', '')
    movie_name_i = movie_name_i.replace('SpotsFiltered', '')
    movie_name_i = movie_name_i.replace('analysis', '')
    movie_name_i = movie_name_i.replace('Analysis', '')
    movie_name_i = movie_name_i.replace('Analyse','')
    movie_name_i = movie_name_i.replace('-Hist-Hetero', '')
    movie_name_i = movie_name_i.replace('-MCP', '')
    movie_name_i = movie_name_i.replace('Hist','')
    movie_name_i = movie_name_i.replace('Homo','')
    movie_name_i = movie_name_i.replace('-','')
    movie_name_i = movie_name_i.replace(' By3.6', '')
    movie_name_i = movie_name_i.replace('Mother','')
    movie_name_i = movie_name_i.replace('Hetero','')
    movie_name_i = movie_name_i.replace('xMCPhomo', '')
    movie_name_i = movie_name_i.replace('xMCP', '')
    movie_name_i = movie_name_i.replace('xMCP', '')
    movie_name_i = movie_name_i.replace('sna24x', '')
    movie_name_i = movie_name_i.replace('noswhite', '')
    movie_name_i = movie_name_i.replace('CF', '')
    movie_name_i = movie_name_i.replace(' S7V5 ','')
    movie_name_i = movie_name_i.replace('IM', '')
    movie_name_i = movie_name_i.replace('Trace', '')      
    movie_name_i = movie_name_i.replace('By7.72s', '')  
    movie_name_i = movie_name_i.replace('by7.72', '') 
    movie_name_i = movie_name_i.replace('  ',' ')
    movie_name_i = movie_name_i.replace('3.65','')
    movie_name_i = movie_name_i.replace('by','')
    movie_name_i = movie_name_i.replace('XMCP','')
    movie_name_i = movie_name_i.replace('MCP','')
    if movie_name_i[-1] == ' ':
        movie_name_i = movie_name_i[:-1]
    return movie_name_i
