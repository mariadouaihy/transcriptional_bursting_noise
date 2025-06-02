#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:14:02 2024

@author: mdouaihy

type_of_analysis = {pospred_per_nuclei, pospred_in_time, ms2_signal_in_time, ms2_signal_per_nuclei}

parameters with star are optional

- pospred_per_nuclei:
    Objective   : here we are calculating the noise coming from the varying of the sum of pol II in differennuclei
                so how long does the noise of a nuclei vary amoung itself across time
    
    Input       :
    1) inputpath    : path to the data
    2) extension    : raw data extension type
    3) phenotype    : the different lines we want to compare
    4) outputpath   : where we want to store the results
    6) fParam       : parameter full path
    7) tfinal       : when we want to finish the analysis so that we have the same length of movie for all nuclei
    8) *data_type   : same as in burstdeconv it specifies if we names the deconvolved files an additional name
    
    
    
- pospred_in_time:
    Objective: noise per time point where for each time point we calculate the number of pol II contributing to this specific time point w.r.t. 
                the dwell time
    
    Input:
    1) inputpath    : path to the data
    2) extension    : raw data extension type
    3) phenotype    : the different lines we want to compare
    4) outputpath   : where we want to store the results
    5) TW_min       : in minutes. Time window between each time points that we are calculating the number of pol II that contributed to this specific time point
    6) fParam       : parameter full path
    7) tfinal       : when we want to finish the analysis so that we have the same length of movie for all nuclei
    8) *data_type   : same as in burstdeconv it specifies if we names the deconvolved files an additional name
    
    
- ms2_signal_in_time:
    Objective:  we take the average in a time window and then calculate the noise across different nuclei/movies for the same time window
    
    Input:
    1) inputpath    : path to the data
    2) extension    : raw data extension type
    3) phenotype    : the different lines we want to compare
    4) outputpath   : where we want to store the results
    5) data_used    : type of data used: DataExp for raw data and DataPred for reconstructed data
    6) TW_min       : time window in minutes
    7) *data_type   : same as in burstdeconv it specifies if we names the deconvolved files an additional name
    8) *fParam      : if using datapred must be added. parameter full path
    9) *calibration : if data is not calibrated add this if you want to add the calibration else if the data is already calibrated then calibration = 1

    
- ms2_signal_per_nuclei:
    Objective: we calculate the noise across different nuclei
    
    Input:
    1) inputpath    : path to the data
    2) extension    : raw data extension type
    3) phenotype    : the different lines we want to compare
    4) outputpath   : where we want to store the results
    5) data_used    : type of data used: DataExp for raw data and DataPred for reconstructed data
    6) *data_type   : same as in burstdeconv it specifies if we names the deconvolved files an additional name
    7) *fParam      : if using datapred must be added. parameter full path
    8) *calibration : if data is not calibrated add this if you want to add the calibration else if the data is already calibrated then calibration = 1

         
"""


import numpy as np
import sys
import os

import matplotlib.pyplot as plt
plt.close('all')

sys.path.append('./utilities/' ) 
from Pre_Noise_Calculations_classes_130524 import PosPredInNuclei, PosPredPreparation, PosPredInTime, TranscriptionalSignalInTime, TranscriptionSignalPreparation, TranscriptionalSignalInNuclei
from movies_combining import movies_combining_rawData



class noise():
    def __init__(self,inputpath, phenotype_gene, output_folder, tfinal=1000, TW_min=1, calibration=1):
        phenotype_rawData       = [i + 'rawData/' for i in phenotype_gene]
        phenotype_BurstDeconv   = [i + 'DeconvOutput/' for i in phenotype_gene]
        fParam_name             = 'HyperParameters.npz'

        """ additional parameters """
        data_used               = 'DataExp'
        extension               = '.xlsx'
        data_type               = ''


        tfinal = self.check_tfinal(tfinal, inputpath, phenotype_rawData, data_type, extension, calibration)



        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        self.pospred_in_nuclei(inputpath, extension, phenotype_gene, phenotype_rawData, phenotype_BurstDeconv, output_folder, fParam_name, tfinal, data_type)
        self.pospred_in_time(inputpath, extension, phenotype_gene, phenotype_rawData, phenotype_BurstDeconv, output_folder, TW_min, fParam_name, tfinal, data_type)
        self.ms2_signal_in_time(inputpath, extension, phenotype_rawData, phenotype_gene, output_folder, data_used, tfinal, TW_min, data_type, phenotype_BurstDeconv, fParam_name, calibration, color_heatmap='YlOrBr')
        self.ms2_signal_per_nuclei(inputpath, extension, phenotype_rawData, phenotype_gene, output_folder, data_used, tfinal, data_type, phenotype_BurstDeconv, fParam_name, calibration, color_heatmap='YlOrBr')



    def check_tfinal(self, tfinal, inputpath, phenotype_rawData, data_type, extension, calibration):

        for gene_i in range(len(phenotype_rawData)):
            inputpath_i                                              = inputpath + phenotype_rawData[gene_i]
            xlsPath                                                  = inputpath_i

            [FrameLen, data_mitosis_combined, tmax_combined, tstart] = movies_combining_rawData(xlsPath, calibration, extension)
            tfinal_i = data_mitosis_combined.shape[0] * FrameLen / 60

            if tfinal_i < tfinal:
                tfinal = tfinal_i

                print('!!! tfinal ' + str(tfinal))

        return tfinal

    def pospred_in_nuclei(self, inputpath, extension, genotype_gene, phenotype_rawData, phenotype_BurstDeconv, output_folder, fParam_name, tfinal, data_type):
        pospred_preparation = PosPredPreparation(inputpath, extension, data_type, genotype_gene, phenotype_rawData, phenotype_BurstDeconv, fParam_name, tfinal)
        """ variation of pol II numbers in a window between nuclei."""

        data_per_pheno_mitosis = pospred_preparation.data_per_pheno_mitosis
        data_per_movie_mitosis = pospred_preparation.data_per_movie_mitosis

        lgd_name_list = pospred_preparation.lgd_name_list
        lgd_per_movie_list = pospred_preparation.lgd_per_movie_list

        FrameLen_list = pospred_preparation.FrameLen_list

        pol_II_nbr_TW_per_movie = []
        pol_II_nbr_TW_per_pheno = []
        for pheno_i in range(len(lgd_name_list)):
            pol_II_nbr_TW_per_pheno.append(np.nansum(data_per_pheno_mitosis[pheno_i], axis=0))

            pol_II_nbr_TW_per_movie_i = []
            for movie_i in data_per_movie_mitosis[pheno_i]:
              pol_II_nbr_TW_per_movie_i.append(np.nansum(movie_i, axis=0))
            pol_II_nbr_TW_per_movie.append(pol_II_nbr_TW_per_movie_i)

        noise_per_nuclei_analys = PosPredInNuclei(pol_II_nbr_TW_per_pheno, pol_II_nbr_TW_per_movie, output_folder, inputpath, 'PosPred' , lgd_name_list, genotype_gene, lgd_per_movie_list, FrameLen_list)


        """movie per movie"""
        for i in range(len(lgd_name_list)):
            df_data  =  noise_per_nuclei_analys.Organize_Noise_per_Embryo_Pol_II(pol_II_nbr_TW_per_movie[i],pol_II_nbr_TW_per_pheno[i], lgd_per_movie_list[i], [FrameLen_list[i]]*len( lgd_per_movie_list[i]), genotype_gene[i])
            noise_per_nuclei_analys.Plot_NoisePerNuclei_Pol_II(df_data, lgd_name_list[i], genotype_gene[i])


        """ all reslts """
        noise_per_nuclei_analys.noise_pospres_total(pol_II_nbr_TW_per_movie, genotype_gene, lgd_name_list)
        xlsfile_name = noise_per_nuclei_analys.output_folder + noise_per_nuclei_analys.outputPath_extension + '/' + noise_per_nuclei_analys.data_used + '_total_noise_Info_per_nuclei.xlsx'
        noise_per_nuclei_analys.plotting_results_total_violonPlot(xlsfile_name)

    ##############################

    def pospred_in_time(self, inputpath, extension, genotype_gene, phenotype_rawData, phenotype_BurstDeconv, output_folder, TW_min, fParam_name, tfinal, data_type):
        pospred_preparation = PosPredPreparation(inputpath, extension, data_type, genotype_gene, phenotype_rawData, phenotype_BurstDeconv, fParam_name, tfinal, TW_min)

        """ extracting pol 2 contributing in a time window. w.r.t. dwell time in time """

        data_per_pheno_mitosis = pospred_preparation.data_per_pheno_mitosis
        data_per_movie_mitosis = pospred_preparation.data_per_movie_mitosis

        time_points_list = pospred_preparation.time_points_list
        lgd_per_movie_list = pospred_preparation.lgd_per_movie_list
        lgd_name_list = pospred_preparation.lgd_name_list





        pospred_in_time_results = PosPredInTime(TW_min, output_folder, inputpath, genotype_gene)


        pol_II_per_TW_movie_mitosis = []
        pol_II_per_TW_pheno_mitosis = []
        pol_II_per_TW_movie_aligned = []


        for pheno_i in range(len(phenotype_BurstDeconv)):
            fParam             = inputpath + phenotype_BurstDeconv[pheno_i] + fParam_name
            parameter_content       = np.load(fParam)

            Polym_speed        = parameter_content['Polym_speed']
            TaillePostMarq     = parameter_content['TaillePostMarq']
            TailleSeqMarq      = parameter_content['TailleSeqMarq']
            FreqEchSimu        = parameter_content['FreqEchSimu']
            signal_length      = int((TailleSeqMarq+TaillePostMarq)/Polym_speed)

            pol_II_per_TW_pheno_mitosis.append(pospred_in_time_results.pol_II_nbr_at_t_wrt_mitosis(data_per_pheno_mitosis[pheno_i], time_points_list[pheno_i]*60-signal_length, TailleSeqMarq, TaillePostMarq, Polym_speed, FreqEchSimu, signal_length))

            pol_II_per_TW_movie_mitosis_i = []
            pol_II_per_TW_movie_aligned_i = []
            for movie_i in range(len(lgd_per_movie_list[pheno_i])):
                pol_II_per_TW_movie_mitosis_i.append(pospred_in_time_results.pol_II_nbr_at_t_wrt_mitosis(data_per_movie_mitosis[pheno_i][movie_i], time_points_list[pheno_i]*60-signal_length, TailleSeqMarq, TaillePostMarq, Polym_speed, FreqEchSimu, signal_length))
                pol_II_per_TW_movie_aligned_i.append(pospred_in_time_results.pol_II_nbr_at_t_wrt_aligned(data_per_movie_mitosis[pheno_i][movie_i], time_points_list[pheno_i]*60-signal_length, TailleSeqMarq, TaillePostMarq, Polym_speed, FreqEchSimu, signal_length))

            pol_II_per_TW_movie_mitosis.append(pol_II_per_TW_movie_mitosis_i)
            pol_II_per_TW_movie_aligned.append(pol_II_per_TW_movie_aligned_i)





        """ saving to data frame and xlsx sheets"""
        pospred_in_time_results.saving_pos_pred_noise(pol_II_per_TW_movie_mitosis, time_points_list, lgd_name_list, lgd_per_movie_list, genotype_gene, pol_II_per_TW_pheno_mitosis)
        pospred_in_time_results.saving_pos_pred_noise(pol_II_per_TW_movie_aligned, time_points_list, lgd_name_list, lgd_per_movie_list, genotype_gene)




        """ plotting """

        lw = 1



        """ Mean w.r.t. averaged error = mean +-std """
        data_averaged_list_mitosis  = [[np.nanmean(i,axis=1) for i in pol_II_per_TW_movie_mitosis[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        data_combined_list_mitosis  = [np.nanmean(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)  for pheno_i in range(len(lgd_name_list))]
        data_averaged_list_aligned  = [[np.nanmean(i,axis=1) for i in pol_II_per_TW_movie_aligned[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        title                       = 'Mean'
        pospred_in_time_results.plotting_pospred_phenotypes(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, title,  lw)
        pospred_in_time_results.plotting_pospred_per_movies(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, genotype_gene, lgd_per_movie_list, title, lw)


        """ Variance w.r.t. averaged error = mean +-std """
        data_averaged_list_mitosis  = [[np.nanvar(i,axis=1) for i in pol_II_per_TW_movie_mitosis[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        data_combined_list_mitosis  = [np.nanvar(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)  for pheno_i in range(len(lgd_name_list))]
        data_averaged_list_aligned  = [[np.nanvar(i,axis=1) for i in pol_II_per_TW_movie_aligned[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        title                       = 'Variance'
        pospred_in_time_results.plotting_pospred_phenotypes(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, title, lw)
        pospred_in_time_results.plotting_pospred_per_movies(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list,  genotype_gene, lgd_per_movie_list, title, lw)


        """ FF w.r.t. averaged error = mean +-std """
        data_averaged_list_mitosis  = [[np.nanvar(i,axis=1)/np.nanmean(i,axis=1) for i in pol_II_per_TW_movie_mitosis[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        data_combined_list_mitosis  = [np.nanvar(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)/np.nanmean(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)  for pheno_i in range(len(lgd_name_list))]
        data_averaged_list_aligned  = [[np.nanvar(i,axis=1)/np.nanmean(i,axis=1) for i in pol_II_per_TW_movie_aligned[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        title                       = 'FF'
        pospred_in_time_results.plotting_pospred_phenotypes(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, title, lw)
        pospred_in_time_results.plotting_pospred_per_movies(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, genotype_gene, lgd_per_movie_list, title, lw)


        """ COV w.r.t. averaged error = mean +-std """
        data_averaged_list_mitosis  = [[np.nanstd(i,axis=1)/np.nanmean(i,axis=1) for i in pol_II_per_TW_movie_mitosis[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        data_combined_list_mitosis  = [np.nanstd(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)/np.nanmean(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)  for pheno_i in range(len(lgd_name_list))]
        data_averaged_list_aligned  = [[np.nanstd(i,axis=1)/np.nanmean(i,axis=1) for i in pol_II_per_TW_movie_aligned[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        title                       = 'COV'
        pospred_in_time_results.plotting_pospred_phenotypes(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, title, lw)
        pospred_in_time_results.plotting_pospred_per_movies(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, genotype_gene, lgd_per_movie_list, title, lw)


        """ CV**2 w.r.t. averaged error = mean +-std """
        data_averaged_list_mitosis  = [[np.nanvar(i,axis=1)/np.nanmean(i,axis=1)**2 for i in pol_II_per_TW_movie_mitosis[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        data_combined_list_mitosis  = [np.nanvar(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)/np.nanmean(pol_II_per_TW_pheno_mitosis[pheno_i], axis=1)**2  for pheno_i in range(len(lgd_name_list))]
        data_averaged_list_aligned  = [[np.nanvar(i,axis=1)/np.nanmean(i,axis=1)**2 for i in pol_II_per_TW_movie_aligned[pheno_i]] for pheno_i in range(len(lgd_name_list))]
        title                       = 'CV**2'
        pospred_in_time_results.plotting_pospred_phenotypes(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, title, lw)
        pospred_in_time_results.plotting_pospred_per_movies(data_averaged_list_mitosis, data_combined_list_mitosis, data_averaged_list_aligned, time_points_list, lgd_name_list, genotype_gene, lgd_per_movie_list, title, lw)

    def ms2_signal_in_time(self, inputpath, extension, phenotype_rawData, genotype_gene, output_folder, data_used, tfinal, TW_min, data_type = '', phenotype_BurstDeconv = '', fParam_name = '', calibration=1, color_heatmap = 'YlOrBr'):

        """ Noise in time windows from transcriptional signal """
        if fParam_name == '':
            ms2_signal_preparation = TranscriptionSignalPreparation(inputpath, data_used, tfinal, TW_min, phenotype_rawData, genotype_gene, data_type, extension, calibration, color_heatmap)
        else:
            ms2_signal_preparation = TranscriptionSignalPreparation(inputpath, data_used, tfinal, TW_min, phenotype_rawData, genotype_gene, data_type, extension, calibration, color_heatmap , phenotype_BurstDeconv, fParam_name)


        data_per_pheno_mitosis = ms2_signal_preparation.data_per_pheno_mitosis
        data_per_movie_mitosis = ms2_signal_preparation.data_per_movie_mitosis
        data_per_pheno_aligned = ms2_signal_preparation.data_per_pheno_aligned
        data_per_movie_aligned = ms2_signal_preparation.data_per_movie_aligned


        FrameLen_list = ms2_signal_preparation.FrameLen_list

        lgd_name_list = ms2_signal_preparation.lgd_name_list
        lgd_per_movie_list = ms2_signal_preparation.lgd_per_movie_list

        transcription_in_time = TranscriptionalSignalInTime(data_per_pheno_mitosis, data_per_movie_mitosis,
                 data_per_pheno_aligned, data_per_movie_aligned,
                 FrameLen_list, lgd_name_list, lgd_per_movie_list,
                 output_folder, inputpath, data_used, TW_min, genotype_gene, color_heatmap)

        """ plotting """
        transcription_in_time.plot_act_per_phenotype(data_per_pheno_mitosis, data_per_movie_mitosis, FrameLen_list,lgd_name_list ) # activation per phenotype
        transcription_in_time.plot_mean(data_per_pheno_mitosis, data_per_pheno_aligned, FrameLen_list, lgd_name_list  ) # mean from combined movies from phenotype
        transcription_in_time.plot_per_phenotype(data_per_pheno_mitosis, data_per_pheno_aligned, data_per_movie_mitosis, data_per_movie_aligned, FrameLen_list, lgd_name_list, 'Variance')
        transcription_in_time.plot_per_phenotype(data_per_pheno_mitosis, data_per_pheno_aligned, data_per_movie_mitosis, data_per_movie_aligned, FrameLen_list, lgd_name_list, 'FF')
        transcription_in_time.plot_per_phenotype(data_per_pheno_mitosis, data_per_pheno_aligned, data_per_movie_mitosis, data_per_movie_aligned, FrameLen_list, lgd_name_list, 'COV')
        transcription_in_time.plot_per_phenotype(data_per_pheno_mitosis, data_per_pheno_aligned, data_per_movie_mitosis, data_per_movie_aligned, FrameLen_list, lgd_name_list, 'CV**2')

        for pheno_i in range(len(data_per_pheno_mitosis)):
            transcription_in_time.plot_act_per_movie(data_per_movie_mitosis[pheno_i], FrameLen_list[pheno_i]*np.ones(len(lgd_per_movie_list[pheno_i])), lgd_per_movie_list[pheno_i] , lgd_name_list[pheno_i], genotype_gene[pheno_i])
            transcription_in_time.plot_mean(data_per_movie_mitosis[pheno_i], data_per_movie_aligned[pheno_i], FrameLen_list[pheno_i]*np.ones(len(lgd_per_movie_list[pheno_i])), lgd_per_movie_list[pheno_i], lgd_name_list[pheno_i], genotype_gene[pheno_i] )
            transcription_in_time.plot_per_movie(data_per_pheno_mitosis[pheno_i], data_per_pheno_aligned[pheno_i], data_per_movie_mitosis[pheno_i], data_per_movie_aligned[pheno_i], FrameLen_list[pheno_i], lgd_per_movie_list[pheno_i], lgd_name_list[pheno_i], 'Variance', genotype_gene[pheno_i])
            transcription_in_time.plot_per_movie(data_per_pheno_mitosis[pheno_i], data_per_pheno_aligned[pheno_i], data_per_movie_mitosis[pheno_i], data_per_movie_aligned[pheno_i], FrameLen_list[pheno_i], lgd_per_movie_list[pheno_i], lgd_name_list[pheno_i], 'FF', genotype_gene[pheno_i])
            transcription_in_time.plot_per_movie(data_per_pheno_mitosis[pheno_i], data_per_pheno_aligned[pheno_i], data_per_movie_mitosis[pheno_i], data_per_movie_aligned[pheno_i], FrameLen_list[pheno_i], lgd_per_movie_list[pheno_i], lgd_name_list[pheno_i], 'COV', genotype_gene[pheno_i])
            transcription_in_time.plot_per_movie(data_per_pheno_mitosis[pheno_i], data_per_pheno_aligned[pheno_i], data_per_movie_mitosis[pheno_i], data_per_movie_aligned[pheno_i], FrameLen_list[pheno_i], lgd_per_movie_list[pheno_i], lgd_name_list[pheno_i], 'CV**2', genotype_gene[pheno_i])

        """ saving results """
        transcription_in_time.save_results_aligned(data_per_pheno_aligned, data_per_movie_aligned, FrameLen_list, lgd_name_list, lgd_per_movie_list, genotype_gene)
        transcription_in_time.save_results_into_motisis(data_per_pheno_mitosis, data_per_movie_mitosis, FrameLen_list, lgd_name_list, lgd_per_movie_list, genotype_gene)

    def ms2_signal_per_nuclei(self, inputpath, extension, phenotype_rawData, genotype_gene, output_folder, data_used, tfinal, data_type = '', phenotype_BurstDeconv = '', fParam_name = '', calibration=1, color_heatmap = 'YlOrBr'):
        TW_min  = 1
        if fParam_name == '':
            ms2_signal_preparation = TranscriptionSignalPreparation(inputpath, data_used, tfinal, TW_min, phenotype_rawData, genotype_gene, data_type, extension, calibration, color_heatmap)
        else:
            ms2_signal_preparation = TranscriptionSignalPreparation(inputpath, data_used, tfinal, TW_min, phenotype_rawData, genotype_gene, data_type, extension, calibration, color_heatmap , phenotype_BurstDeconv, fParam_name)


        data_per_pheno_mitosis = ms2_signal_preparation.data_per_pheno_mitosis
        data_per_movie_mitosis = ms2_signal_preparation.data_per_movie_mitosis

        FrameLen_list = ms2_signal_preparation.FrameLen_list

        lgd_name_list = ms2_signal_preparation.lgd_name_list
        lgd_per_movie_list = ms2_signal_preparation.lgd_per_movie_list

        transcription_in_nuclei = TranscriptionalSignalInNuclei(data_per_pheno_mitosis, data_per_movie_mitosis, output_folder, inputpath,  data_used, lgd_name_list, genotype_gene, lgd_per_movie_list, FrameLen_list)

        """ Noise per nuclei from transcriptional signal """
        for i in range(len(lgd_name_list)):
            df_data = transcription_in_nuclei.Organize_Noise_per_Nuclei(data_per_movie_mitosis[i],lgd_per_movie_list[i], [FrameLen_list[i]]*len( lgd_per_movie_list[i]))
            transcription_in_nuclei.Save_NoisePernuclei(df_data, data_per_pheno_mitosis[i],  genotype_gene[i], FrameLen_list[i])
            transcription_in_nuclei.Plot_NoisePerNuclei(df_data, lgd_name_list[i], genotype_gene[i])

        transcription_in_nuclei.noise_datapred_total(data_per_movie_mitosis, genotype_gene, lgd_name_list)

        xlsfile_name = transcription_in_nuclei.output_folder + transcription_in_nuclei.outputPath_extension + '/' + transcription_in_nuclei.data_used + '_total_noise_Info_per_nuclei.xlsx'
        transcription_in_nuclei.plotting_results_total_violonPlot(xlsfile_name)

        xlsfile_name = transcription_in_nuclei.output_folder + transcription_in_nuclei.outputPath_extension + '/AUC_' + transcription_in_nuclei.data_used + '_total_noise_Info_per_nuclei.xlsx'
        transcription_in_nuclei.plotting_results_total_violonPlot(xlsfile_name, 'AUC_')

        plt.close('all')

noise("./data/", ["twi-mDPE/", "twi-mDPEpTATA/"], "./noise_output/")