#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 10:27:39 2024

@author: mdouaihy

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys
from matplotlib.patches import Patch

sys.path.append('./utilities/')

from movies_combining import movies_combining, data_reshape, clean_movie_name, data_average_in_TW, movies_combining_rawData, allign_data_Activation, data_reshape_raw
from extentForPlot import extentForPlot

plt.close('all')
class TranscriptionalSignalInTime():
    def __init__(self, data_per_pheno_mitosis, data_per_movie_mitosis,
                 data_per_pheno_aligned, data_per_movie_aligned,
                 FrameLen_list, lgd_name_list, lgd_per_movie_list,
                 output_folder, inputpath, data_used, TW_min, phenotype_gene, color_heatmap = 'YlOrBr'):
        
        self.output_folder             = output_folder
        self.outputPath_extension   = data_used + '_by_time_TW_' + str(TW_min) + '/'
        self.inputpath              = inputpath
        self.TW_min                 = TW_min
        self.color_heatmap          = color_heatmap
        
        self.data_per_pheno_mitosis = data_per_pheno_mitosis
        self.data_per_movie_mitosis = data_per_movie_mitosis
        
        self.data_per_pheno_aligned = data_per_pheno_aligned
        self.data_per_movie_aligned = data_per_movie_aligned
        
        self.FrameLen_list          = FrameLen_list    
        self.lgd_name_list          = lgd_name_list
        self.lgd_per_movie_list     = lgd_per_movie_list
        
        self.data_used              = data_used
        
        """ create outputFolder """
        outputPath = self.output_folder + self.outputPath_extension
        if not os.path.exists(outputPath): 
            os.mkdir(outputPath)
            
        

    """ percentage of active nuclei. Only from mitosis """
    def percentage_act(self, data):
        return np.count_nonzero(data, axis=1) / np.count_nonzero(~np.isnan(data), axis=1) * 100
    

            

        
    def plot_act_per_phenotype(self, data_per_pheno, data_per_movie, FrameLen_list, lgd_name_list):
        h, ax           = plt.subplots()
        resolution = 1.5
        for gene_i in range(len(data_per_pheno)):
            color = plt.cm.jet(gene_i / ((len(data_per_pheno)+1) * resolution))  # Adjust colormap as needed
            TW_bp               = int(self.TW_min*60/FrameLen_list[gene_i])
            data                = data_average_in_TW(data_per_pheno[gene_i], TW_bp)
            X_TW_min            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen_list[gene_i]/60
            ax.plot(X_TW_min,  self.percentage_act(data), label = lgd_name_list[gene_i], lw = 3, color = color)
            
            data                = [data_average_in_TW(i, TW_bp) for i in data_per_movie[gene_i]]
            error               = np.array([self.percentage_act(j[:len(X_TW_min)]) for j in data])
            ax.fill_between(X_TW_min, np.nanmin(error, axis=0), np.nanmax(error, axis=0), alpha = 0.2, color = color) 
        
        ax.legend()
        ax.set_xlabel('Time [min]', fontsize=12)
        ax.set_title('% non zero nuclei TW = ' + str(self.TW_min) + ' min', fontsize=12)
        h.savefig(self.output_folder + self.outputPath_extension + 'ActiveNuclei_TW_' + str(self.TW_min) + '.pdf')
        
    def plot_act_per_movie(self, data_per_movie, FrameLen_list, lgd_name_list, lgd_pheno, phenotype_path):
        h, ax           = plt.subplots()
        
        resolution = 1.5
        for gene_i in range(len(data_per_movie)):
            color = plt.cm.jet(gene_i / ((len(data_per_movie)+1) * resolution))  # Adjust colormap as needed
            TW_bp               = int(self.TW_min*60/FrameLen_list[gene_i])
            data                = data_average_in_TW(data_per_movie[gene_i], TW_bp)
            X_TW_min            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen_list[gene_i]/60
            ax.plot(X_TW_min,  self.percentage_act(data), label = lgd_name_list[gene_i], lw = 3,  color = color)

        ax.legend()
        ax.set_xlabel('Time [min]', fontsize=12)
        ax.set_title('% non zero nuclei ' + lgd_pheno.replace(' ','').replace('_','') + ': TW = ' + str(self.TW_min) + ' min', fontsize=12)
        h.savefig(self.output_folder + self.outputPath_extension + 'ActiveNuclei_TW_' + str(self.TW_min) +  '_' + lgd_pheno.replace(' ','').replace('_','') + '.pdf')
        
            
    """ mean of data with +- standard deviation """
    
    def plot_mean(self, Data_mitosis_list, Data_aligned_list, FrameLen_list, lgd_name_list,  lgd_pheno = '' , phenotype_path = ''):
        h, ax           = plt.subplots(2,1, figsize = (7,9))
        resolution = 1.5
        for gene_i in range(len(Data_mitosis_list)):
            color = plt.cm.jet(gene_i / ((len(Data_mitosis_list)+1) * resolution))  # Adjust colormap as needed
            TW_bp               = int(self.TW_min*60/FrameLen_list[gene_i])
            data                = data_average_in_TW(Data_mitosis_list[gene_i], TW_bp)
            
            X_TW_min            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen_list[gene_i]/60        
            ax[0].plot(X_TW_min, np.nanmean(data, axis = 1), label = lgd_name_list[gene_i].replace(' ','').replace('_',''), lw=3, color = color)
            ax[0].fill_between(X_TW_min, np.nanmean(data, axis = 1) -  np.nanstd(data, axis = 1),np.nanmean(data, axis = 1) +  np.nanstd(data, axis = 1), alpha = 0.2, color = color) 
            ax[0].legend()
            ax[0].set_xlabel('Time [min]', fontsize=12)
    
    
            data                = data_average_in_TW(Data_aligned_list[gene_i], TW_bp)
            X_TW_min            = np.arange(0, data.shape[0])
            ax[1].plot(X_TW_min, np.nanmean(data, axis = 1), label = lgd_name_list[gene_i].replace(' ','').replace('_',''), lw=3, color = color)
            ax[1].fill_between(X_TW_min, np.nanmean(data, axis = 1) -  np.nanstd(data, axis = 1),np.nanmean(data, axis = 1) +  np.nanstd(data, axis = 1), alpha = 0.2, color = color) 
            ax[1].legend()
            ax[1].set_xlabel('Time [min]', fontsize=12)
            

        ax[0].set_title('From Mitosis')
        ax[1].set_title('aligned w.r.t. 20% of max intensity')
        h.suptitle('Mean : ' +  lgd_pheno.replace(' ','').replace('_','') + ' TW = ' + str(self.TW_min) + ' min', fontsize=12)    
        plt.tight_layout()
        
        h.savefig(self.output_folder + self.outputPath_extension + 'mean_TW_' + str(self.TW_min) + '_' + lgd_pheno.replace(' ','').replace('_','') + '.pdf')

    """ heatmap of dataexp """
    
    def plot_HeatMap_DataExp(self, data, lgd_name, FrameLen, tmax_combined):
        
        X_min           = np.arange(0, data.shape[0])*FrameLen/60
        cm_jet          = plt.cm.get_cmap(self.color_heatmap)
        
        ### compute T0 at 10% of maximum intensity                
        T0_smoothed     = np.zeros((data.shape[1],))  ##### will contain the start of the analyzed region
        for i in range(data.shape[1]): ### for all cells
            max_intensity       = np.nanmax(data)
            ihit                = np.where(data[:,i] > max_intensity/10 )[0]
            if len(ihit)!=0:
                ihi_min                 =  min(ihit)
                T0_smoothed[i]          = (ihi_min-1)*FrameLen; #### T0_combined    
            else:
                T0_smoothed[i]          = tmax_combined[i]
            
    
        
        argsort         = T0_smoothed.argsort()
        data_sorted     = data[:,argsort]
        
        h               = plt.figure()
        Y_normal        = np.arange(1,data.shape[1]+1)
        plt.imshow(np.nan_to_num(data_sorted).T, cmap=cm_jet, extent = extentForPlot(X_min).result + extentForPlot(Y_normal[::-1]).result, aspect='auto', origin='upper')
        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('Transcription site', fontsize=12)
        plt.colorbar()
        plt.title(lgd_name)
        figfile         = self.output_folder + self.outputPath_extension + self.data_used + lgd_name.replace(' ','_') + '.pdf'
        h.savefig(figfile)
        
    
    """ variance of data with +- min max """
    

    def mean_calculation(self, data):
        return np.nanmean(data, axis=1)

    def variance_calculation(self, data):
        return np.nanvar(data, axis=1)

       
    def FF_calculation(self, data):
        return np.nanvar(data, axis=1)/np.nanmean(data, axis=1)

    def COV_calculation(self, data):
        return np.nanstd(data, axis=1)/np.nanmean(data, axis=1)

    def CV_calculation(self, data):
        return np.var(data, axis=1)/np.nanmean(data, axis=1)**2

    
    def plot_per_phenotype(self, data_per_pheno_mitosis, data_per_pheno_aligned, data_per_movie_mitosis, data_per_movie_aligned, FrameLen_list, lgd_name_list, title):
        resolution = 1.5
        
        title_to_function = {'Mean': self.mean_calculation, 'Variance': self.variance_calculation,
                            'FF': self.FF_calculation, 'COV': self.COV_calculation, 'CV**2': self.CV_calculation}

        statistics = title_to_function.get(title)
   


        h, ax           = plt.subplots(2,1, figsize = (7,9))
        
        for gene_i in range(len(data_per_pheno_mitosis)):
            color = plt.cm.jet(gene_i / ((len(data_per_pheno_mitosis)+1) * resolution))  # Adjust colormap as needed
            
            """ mitosis """
            TW_bp               = int(self.TW_min*60/FrameLen_list[gene_i])
            data                = data_average_in_TW(data_per_pheno_mitosis[gene_i], TW_bp)
            X_TW_min            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen_list[gene_i]/60
            
            data_per_movie      = np.array([statistics(data_average_in_TW(j, TW_bp)[:len(X_TW_min),:]) for j in data_per_movie_mitosis[gene_i]])
#            data_combined       = statistics(data)
            data_averaged       = np.nanmean(data_per_movie, axis=0)
            
            ax[0].fill_between(X_TW_min, np.nanmin(data_per_movie, axis=0), np.nanmax(data_per_movie, axis=0), alpha = 0.2, color=color) 
#            ax[0].plot(X_TW_min,  data_combined, label =  lgd_name_list[gene_i] + ' Combined', lw = 3, alpha = 0.4, color=color)
            ax[0].plot(X_TW_min,  data_averaged, label = lgd_name_list[gene_i] + ' Averaged', lw = 3, color=color)
            
            
            """ alligned """
            data                = data_average_in_TW(data_per_pheno_aligned[gene_i], TW_bp)
            X_TW_min            = np.arange(0, data.shape[0])
            
            data_per_movie      = np.zeros((len(data_per_movie_aligned[gene_i]), len(X_TW_min)))
            data_per_movie      = np.full_like(data_per_movie, np.nan)
            for j_i in  range(len(data_per_movie_aligned[gene_i])):
                j                                            = data_average_in_TW(data_per_movie_aligned[gene_i][j_i], TW_bp)
                data_per_movie_j                             = statistics(j[:len(X_TW_min),:])
                data_per_movie[j_i,:len(data_per_movie_j)]   = data_per_movie_j
#            data_combined       = statistics(data) 
            data_averaged       = np.nanmean(data_per_movie, axis=0)
            
            ax[1].fill_between(X_TW_min, np.nanmin(data_per_movie, axis=0), np.nanmax(data_per_movie, axis=0), alpha = 0.2, color=color) 
#            ax[1].plot(X_TW_min,  data_combined, label = lgd_name_list[gene_i]  + ' Combined', lw = 3, alpha = 0.4, color=color)
            ax[1].plot(X_TW_min,  data_averaged, label = lgd_name_list[gene_i]  + ' Averaged', lw = 3, color=color)
            

        ax[0].legend()
        ax[0].set_xlabel('Time [min]', fontsize=12)
        ax[0].set_title('From Mitosis')

        ax[1].legend()
        ax[1].set_xlabel('Time [min]', fontsize=12)
        ax[1].set_title('aligned w.r.t. 20% of max intensity')
        
        h.suptitle(title + ' : TW = ' + str(self.TW_min) + ' min', fontsize=12)
        
        if (title == 'COV') or (title=='CV**2'):
            ax[1].set_yscale('log')
            ax[0].set_yscale('log')
        plt.tight_layout()
        
        h.savefig(self.output_folder + self.outputPath_extension + title + '_TW_' + str(self.TW_min) + '.pdf')
        
        
        
    def plot_per_movie(self, data_per_pheno_mitosis_i, data_per_pheno_aligned_i, data_per_movie_mitosis, data_per_movie_aligned, FrameLen, lgd_name_list, lgd_pheno, title, phenotype_path):
        resolution = 1.5
        title_to_function           = {'Mean': self.mean_calculation, 'Variance': self.variance_calculation,
                                       'FF': self.FF_calculation, 'COV': self.COV_calculation, 'CV**2': self.CV_calculation}
        statistics                  = title_to_function.get(title)
    
        
        h, ax                       = plt.subplots(2,1, figsize = (7,9))
        
        TW_bp                       = int(self.TW_min*60/FrameLen)
        
        data                        = data_average_in_TW(data_per_pheno_mitosis_i, TW_bp)
        X_TW_min_mitosis            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen/60
        data_combined_mitosis       = statistics(data)        
        data_average_list_mitosis   = []

        data                        = data_average_in_TW(data_per_pheno_aligned_i, TW_bp)
        X_TW_min_aligned            = np.arange(0, data.shape[0])  
        data_combined_aligned       = statistics(data)               
        data_average_list_aligned   = []
        
        for gene_i in range(len(data_per_movie_mitosis)):
            color =  plt.cm.jet(gene_i / ((len(data_per_movie_mitosis)+1) * resolution))
            data                = data_average_in_TW(data_per_movie_mitosis[gene_i], TW_bp)
            X_TW_min            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen/60
            
            data_per_movie      = statistics(data)
            data_average_list_mitosis.append(data_per_movie)
            ax[0].plot(X_TW_min,  data_per_movie, label = lgd_name_list[gene_i], lw = 3,  color = color)
            
            data                = data_average_in_TW(data_per_movie_aligned[gene_i], TW_bp)
            X_TW_min            = np.arange(0, data.shape[0])   
            
            data_per_movie      = statistics(data)
            data_average_list_aligned.append(data_per_movie)
            ax[1].plot(X_TW_min,  data_per_movie, label = lgd_name_list[gene_i], lw = 3, color = color)
        
        color = plt.cm.jet((gene_i+1)/((len(data_per_movie_mitosis)+1) * resolution))  # Adjust colormap as needed
        
        data_averaged       = np.nanmean(np.array([i[:len(X_TW_min_mitosis)] for i in data_average_list_mitosis]), axis = 0)
        ax[0].plot(X_TW_min_mitosis,  data_combined_mitosis, label = 'Combined ', lw = 3, alpha = 0.4, color=color)
        ax[0].plot(X_TW_min_mitosis,  data_averaged, label = 'Averaged ', lw = 3, color=color)

        ax[0].legend()
        ax[0].set_xlabel('Time [min]', fontsize=12)
        ax[0].set_title('From Mitosis')
        
        min_data_len        = min([i.shape[0] for i in data_average_list_aligned])
        data_averaged       = np.nanmean(np.array([i[:min_data_len] for i in data_average_list_aligned]), axis = 0)

        ax[1].plot(X_TW_min_aligned,  data_combined_aligned, label = 'Combined ', lw = 3, alpha = 0.4, color=color)
        ax[1].plot(np.arange(min_data_len),  data_averaged, label = 'Averaged ', lw = 3, color=color)    
        ax[1].legend()
        ax[1].set_xlabel('Time [min]', fontsize=12)
        ax[1].set_title('aligned w.r.t. 20% of max intensity')
        
        h.suptitle(title +  ' : ' + lgd_pheno.replace(' ','').replace('_','')  + ' TW = ' + str(self.TW_min) + ' min', fontsize=12)
        if (title == 'COV') or (title=='CV**2'):
            ax[1].set_yscale('log')
            ax[0].set_yscale('log')
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + title + '_TW_' + str(self.TW_min) +  '_' + lgd_pheno.replace(' ','').replace('_','') + '.pdf')
               
    
    
    """ save results into motisis """
    def save_results_into_motisis(self, data_per_pheno_mitosis, data_per_movie_mitosis, FrameLen_list,  lgd_name_list, lgd_per_movie_list, phenotype_gene):
        
        xlsfilename     = self.output_folder + self.outputPath_extension + self.data_used + '_TW_' + str(self.TW_min) + '_mitosis.xlsx'
        writer          = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
        
            
        for gene_i in range(len(phenotype_gene)):
            TW_bp               = int(self.TW_min*60/FrameLen_list[gene_i])
            data                = data_average_in_TW(data_per_pheno_mitosis[gene_i], TW_bp)
            X_TW_min            = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen_list[gene_i]/60
            
            
            xlsfilename_i     = self.output_folder + self.outputPath_extension + self.data_used + "_" +phenotype_gene[gene_i].replace('/', '') + '_TW_' + str(self.TW_min) + '_mitosis.xlsx'
            writer_i          = pd.ExcelWriter(xlsfilename_i, engine='xlsxwriter')
        
            
            df                  = pd.DataFrame({'Time':  ['T = '+ str(np.round(X_TW_min[i],2)) + ' s' for i in range(len(X_TW_min))],
                                                'Mean': np.nanmean(data, axis = 1),
                                                'Variance': np.nanvar(data, axis = 1),
                                                'FF': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1),
                                                'COV': np.nanstd(data, axis = 1)/np.nanmean(data, axis = 1),
                                                'CV**2': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1)**2})
            df.to_excel(writer, (lgd_name_list[gene_i] + ', combined')[:31], index=False, startrow = 0, startcol = 0)  
            df.to_excel(writer_i, (lgd_name_list[gene_i] + ', combined')[:31], index=False, startrow = 0, startcol = 0)  
            
            width_cell_i        = len('Variance  ')
            writer.sheets[(lgd_name_list[gene_i] + ', combined')[:31]].set_column( 1,5, width_cell_i)
            writer_i.sheets[(lgd_name_list[gene_i] + ', combined')[:31]].set_column( 1,5, width_cell_i)

        
            for movie_i in range(len(data_per_movie_mitosis[gene_i])):
                data                    = data_average_in_TW(data_per_movie_mitosis[gene_i][movie_i], TW_bp)
                X_TW_min                = np.arange(TW_bp/2,  data.shape[0]*TW_bp , TW_bp)*FrameLen_list[gene_i]/60
                
                df                      = pd.DataFrame({'Time':  ['T = '+ str(np.round(X_TW_min[i],2)) + ' s' for i in range(len(X_TW_min))],
                                                        'Mean': np.nanmean(data, axis = 1),
                                                        'Variance': np.nanvar(data, axis = 1),
                                                        'FF': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1),
                                                        'COV': np.nanstd(data, axis = 1)/np.nanmean(data, axis = 1),
                                                        'CV**2': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1)**2})
                df.to_excel(writer_i, sheet_name=(lgd_name_list[gene_i] + ' ' + lgd_per_movie_list[gene_i][movie_i])[:31], index=False, startrow = 0, startcol = 0)  
                
                width_cell_i            = len('Variance  ')
                writer_i.sheets[(lgd_name_list[gene_i] + ' ' + lgd_per_movie_list[gene_i][movie_i])[:31]].set_column( 1,5, width_cell_i)        
            writer_i.close()
        writer.close()
    
    
    
    
    
    
    
    """ save results aligned """
    def save_results_aligned(self, data_per_pheno_aligned, data_per_movie_aligned, FrameLen_list,  lgd_name_list, lgd_per_movie_list, phenotype_gene):
        
        xlsfilename     = self.output_folder + self.outputPath_extension + self.data_used +  '_TW_' + str(self.TW_min) + '_aligned.xlsx'
        writer          = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
        
            
        for gene_i in range(len(phenotype_gene)):
            xlsfilename_i     = self.output_folder + self.outputPath_extension + self.data_used + "_" + phenotype_gene[gene_i].replace('/', '') + '_TW_' + str(self.TW_min) + '_aligned.xlsx'
            writer_i          = pd.ExcelWriter(xlsfilename_i, engine='xlsxwriter')
        
        
        
            TW_bp               = int(self.TW_min*60/FrameLen_list[gene_i])
            data                = data_average_in_TW(data_per_pheno_aligned[gene_i], TW_bp)
            X_TW_min            = np.arange(0, data.shape[0]) 
            
            df                  = pd.DataFrame({'Time':  ['T = '+ str(np.round(X_TW_min[i],2)) + ' s' for i in range(len(X_TW_min))],
                                                'Mean': np.nanmean(data, axis = 1),
                                                'Variance': np.nanvar(data, axis = 1),
                                                'FF': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1),
                                                'COV': np.nanstd(data, axis = 1)/np.nanmean(data, axis = 1),
                                                'CV**2': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1)**2})
            df.to_excel(writer, sheet_name=(lgd_name_list[gene_i] + ', combined')[:31], index=False, startrow = 0, startcol = 0)  
            df.to_excel(writer_i, sheet_name=(lgd_name_list[gene_i] + ', combined')[:31], index=False, startrow = 0, startcol = 0)  
            
            width_cell_i        = len('Variance  ')
            writer.sheets[(lgd_name_list[gene_i] + ', combined')[:31]].set_column( 1,5, width_cell_i)
            writer_i.sheets[(lgd_name_list[gene_i] + ', combined')[:31]].set_column( 1,5, width_cell_i)
            
            for movie_i in range(len(data_per_movie_aligned[gene_i])):
                data                    = data_average_in_TW(data_per_movie_aligned[gene_i][movie_i], TW_bp)
                X_TW_min                = np.arange(0, data.shape[0])
                
                
                df                      = pd.DataFrame({'Time':  ['T = '+ str(np.round(X_TW_min[i],2)) + ' s' for i in range(len(X_TW_min))],
                                                        'Mean': np.nanmean(data, axis = 1),
                                                        'Variance': np.nanvar(data, axis = 1),
                                                        'FF': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1),
                                                        'COV': np.nanstd(data, axis = 1)/np.nanmean(data, axis = 1),
                                                        'CV**2': np.nanvar(data, axis = 1)/np.nanmean(data, axis = 1)**2})
                df.to_excel(writer_i, sheet_name=(lgd_name_list[gene_i] + ' ' + lgd_per_movie_list[gene_i][movie_i])[:31], index=False, startrow = 0, startcol = 0)  
                
                width_cell_i            = len('Variance  ')
                writer_i.sheets[(lgd_name_list[gene_i] + ' ' + lgd_per_movie_list[gene_i][movie_i])[:31]].set_column( 1,5, width_cell_i)        
            writer_i.close()
        writer.close()
    
    
########################################################################################################
















class PosPredInTime():
    def __init__(self, TW_min, output_folder, inputpath, phenotype_gene):
        
        self.TW_min                  = TW_min
        self.output_folder           = output_folder
        self.inputpath               = inputpath
        self.outputPath_extension    = 'PosPred_by_time_TW_' + str(self.TW_min) + '/'
        
        outputPath                   = self.output_folder + self.outputPath_extension
        if not os.path.exists(outputPath): 
            os.mkdir(outputPath)
            

            
    def pol_II_nbr_at_t_wrt_mitosis(self, PosPred, time_points_selected, TailleSeqMarq, TaillePostMarq, Polym_speed, FreqEchSimu, signal_length):

        nn                          = PosPred.shape
        pol_2_nbr_gene_i            = np.zeros((len(time_points_selected),nn[1]))

        for nuclei_i in range(nn[1]):

            pos_nuclei_i                    = np.where(PosPred[:,nuclei_i] == 1)[0]
            times                           = ( pos_nuclei_i +1) / FreqEchSimu - (TailleSeqMarq+TaillePostMarq)/Polym_speed
            for iii in range(len(time_points_selected)):
                time_pt                                 = time_points_selected[iii]
                ind_pol_in_interval                     = np.where((time_pt-signal_length<times) & (times<time_pt))[0] # index of pol II contributed to this point
                pol_2_nbr_gene_i[iii, nuclei_i]         = len(ind_pol_in_interval)
    
        return pol_2_nbr_gene_i
    
    
    def pol_II_nbr_at_t_wrt_aligned(self, PosPred, time_points_selected,
                                    TailleSeqMarq, TaillePostMarq, Polym_speed, FreqEchSimu, signal_length):
        nn                          = PosPred.shape   
        pol_2_nbr_gene_i            = np.zeros((len(time_points_selected),nn[1]))        
        for nuclei_i in range(nn[1]):
            
            pos_nuclei_i                    = np.where(PosPred[:,nuclei_i]==1)[0]
            """ here is the part where we take pospred w.r.t first time ==1 """
            if sum(pos_nuclei_i) ==0:
                PosPred_modified                        = PosPred
            else:
                PosPred_modified                        = PosPred[pos_nuclei_i[0]:,nuclei_i]
                
            ############################
            pos_nuclei_i_mod        = np.where(PosPred_modified == 1)[0]
            times                   = ( pos_nuclei_i_mod +1) / FreqEchSimu -(TailleSeqMarq+TaillePostMarq)/Polym_speed
            for iii in range(len(time_points_selected)):
                time_pt                                 = time_points_selected[iii]
                ind_pol_in_interval                     = np.where((time_pt-signal_length<times) & (times<time_pt))[0]
                pol_2_nbr_gene_i[iii, nuclei_i]         = len(ind_pol_in_interval)
                
        return pol_2_nbr_gene_i
    
    
    
    def plotting_pospred_per_movies(self, data_averaged_list_mitosis, data_combined_list_mitosis, 
                                    data_averaged_list_aligned,
                                    time_points_list, lgd_name_list, phenotype_gene, lgd_per_movie_list, title, lw):
        resolution = 1.5
        for pheno_i in range(len(lgd_name_list)):
            
            h                               = plt.figure(figsize=(7.27, 9.69))
            ax                              = h.subplots(2,1)
        
            times                           = time_points_list[pheno_i]
            
            # averaged
            data                            = data_averaged_list_mitosis[pheno_i]        
            
            for movie_i in range(len(lgd_per_movie_list[pheno_i])):
                color = plt.cm.jet((movie_i)/((len(lgd_per_movie_list[pheno_i])+1) * resolution))  #
                ax[0].plot(times, data[movie_i], label = lgd_per_movie_list[pheno_i][movie_i],  linewidth = 3, color = color)
          
            color = plt.cm.jet((movie_i+1)/((len(lgd_per_movie_list[pheno_i])+1) * resolution))  #
            ax[0].plot(times, np.nanmean(data, axis=0), color = color, label = 'Averaged', linewidth = 3)
            
            # combined
            data                            = data_combined_list_mitosis[pheno_i]
#            ax[0].plot(times, data, color = color, label = 'Combined', linewidth = 3, alpha = 0.4)
            ax[0].set_xlabel('time [min]')
            ax[0].legend()
            
    
            # averaged
            data                            = data_averaged_list_aligned[pheno_i]        
            
            for movie_i in range(len(lgd_per_movie_list[pheno_i])):
                color = plt.cm.jet((movie_i)/((len(lgd_per_movie_list[pheno_i])+1) * resolution))  #
                ax[1].plot(times, data[movie_i], label = lgd_per_movie_list[pheno_i][movie_i], linewidth = 3, color = color )
            
            color = plt.cm.jet((movie_i+1)/((len(lgd_per_movie_list[pheno_i])+1) * resolution))  #
            ax[1].plot(times, np.nanmean(data, axis=0), color = color, label = 'Averaged', linewidth = 3)
            ax[1].legend()
            ax[0].set_title('Time into Mitosis')
            ax[1].set_title('Time into 1st activation')
            h.suptitle(title + ' : ' + lgd_name_list[pheno_i]) 
            if (title == 'COV') or (title=='CV**2'):
                ax[1].set_yscale('log')
                ax[0].set_yscale('log')            
            plt.tight_layout()

            
            h.savefig(self.output_folder + self.outputPath_extension + title + "_" + lgd_name_list[pheno_i] + '.pdf')
            
            
    def plotting_pospred_phenotypes(self, data_averaged_list_mitosis, data_combined_list_mitosis, 
                                    data_averaged_list_aligned,
                                    time_points_list, lgd_name_list, title, lw):
        
        resolution = 1.5
        
        h                           = plt.figure(figsize=(7.27, 9.69))
        ax                          = h.subplots(2,1)
        
        for pheno_i in range(len(lgd_name_list)):
            times                           = time_points_list[pheno_i]
            color = plt.cm.jet(pheno_i / ((len(lgd_name_list)+1) * resolution)) 
            # averaged
            data                            = data_averaged_list_mitosis[pheno_i]     
            ax[0].plot(times, np.nanmean(data, axis=0), color =color, label = lgd_name_list[pheno_i]+ ' Averaged', linewidth = 3 )
            ax[0].fill_between(times, np.nanmean(data, axis=0)-np.nanstd(data, axis=0), np.nanmean(data, axis=0) + np.nanstd(data, axis=0), color = color, alpha = 0.2)
            
            # combined
            data                            = data_combined_list_mitosis[pheno_i]
#            ax[0].plot(times, data, color =color, label = lgd_name_list[pheno_i] + ' combined', linewidth = 3, alpha = 0.4 )
            
            
    
            # averaged
            data                            = data_averaged_list_aligned[pheno_i]        
            ax[1].plot(times, np.nanmean(data, axis=0), color =color, label = lgd_name_list[pheno_i], linewidth = 3 )
            ax[1].fill_between(times, np.nanmean(data, axis=0)-np.nanstd(data, axis=0), np.nanmean(data, axis=0) + np.nanstd(data, axis=0), color = color, alpha = 0.2)
            
    
        ax[0].set_xlabel('time [min]')
        ax[0].legend()
        ax[1].legend()
        
        ax[0].set_title('Time into Mitosis')
        ax[1].set_title('Time into 1st activation')    
        
        h.suptitle(title) 
        
        if (title == 'COV') or (title=='CV**2'):
            ax[1].set_yscale('log')
            ax[0].set_yscale('log')
        plt.tight_layout()
        h.savefig(self.output_folder +  self.outputPath_extension + title + '.pdf')
        
    
    def saving_pos_pred_noise(self, data_per_movie, time_points_list, lgd_name_list, lgd_per_movie_list, phenotype_gene, data_per_pheno = ''):

        if len(data_per_pheno) !=0:
            xlsfilename_all                      = self.output_folder  + self.outputPath_extension + '/PosPred_TW_' + str(self.TW_min) + '_mitosis.xlsx'
            
        else:
            xlsfilename_all                      = self.output_folder  + self.outputPath_extension + '/PosPred_TW_' + str(self.TW_min) + '_aligned.xlsx'
        writer_all = pd.ExcelWriter(xlsfilename_all, engine='xlsxwriter')
     



        for pheno_i in range(len(lgd_name_list)):
            if len(data_per_pheno) !=0:
                xlsfilename                     = self.output_folder + self.outputPath_extension + '/PosPred_' + lgd_name_list[pheno_i] + '_TW_' + str(self.TW_min) + '_mitosis.xlsx'
                
            else:
                xlsfilename                     = self.output_folder + self.outputPath_extension + '/PosPred_' + lgd_name_list[pheno_i] + '_TW_' + str(self.TW_min) + '_aligned.xlsx'
            writer = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
        
            if len(data_per_pheno) !=0:
                """ combined """
                df                              = pd.DataFrame({'Time'       : ['T = '+ str(i) + ' min' for i in time_points_list[pheno_i]],
                                                                'Mean'       : np.nanmean(data_per_pheno[pheno_i], axis=1),
                                                                'Variance'   : np.nanvar(data_per_pheno[pheno_i], axis=1),
                                                                'FF'         : np.nanvar(data_per_pheno[pheno_i], axis=1)/np.nanmean(data_per_pheno[pheno_i], axis=1),
                                                                'COV'        : np.nanstd(data_per_pheno[pheno_i], axis=1)/np.nanmean(data_per_pheno[pheno_i], axis=1),
                                                                'CV**2'      : np.nanvar(data_per_pheno[pheno_i], axis=1)/np.nanmean(data_per_pheno[pheno_i], axis=1)**2,
                                                           })
                df.to_excel(writer, sheet_name= 'Combined', index=False, startrow = 0, startcol = 0)
        
                df.to_excel(writer_all, sheet_name= 'Combined ' + lgd_name_list[pheno_i], index=False, startrow = 0, startcol = 0)
            
            """ Averaged """
            data_mean_per_film              = np.array([np.nanmean(i,axis=1) for i in data_per_movie[pheno_i]])
            data_var_per_film               = np.array([np.nanvar(i,axis=1) for i in data_per_movie[pheno_i]])
            
            df                              = pd.DataFrame({'Time'       :  ['T = '+ str(i) + ' min' for i in time_points_list[pheno_i]],
                                                            'Mean'       : np.nanmean(data_mean_per_film, axis=0),
                                                            'Mean LB'    : np.nanmean(data_mean_per_film, axis=0) - np.nanstd(data_mean_per_film, axis=0),
                                                            'Mean UB'    : np.nanmean(data_mean_per_film, axis=0) + np.nanstd(data_mean_per_film, axis=0),
                               
                                                            'Variance'   : np.nanmean(data_var_per_film, axis=0),
                                                            'Variance LB': np.nanmean(data_var_per_film, axis=0) - np.nanstd(data_var_per_film, axis=0),
                                                            'Variance UB': np.nanmean(data_var_per_film, axis=0) + np.nanstd(data_var_per_film, axis=0),
        
                                                            'FF'         : np.nanmean(data_var_per_film/data_mean_per_film, axis=0),
                                                            'FF LB'      : np.nanmean(data_var_per_film/data_mean_per_film, axis=0) - np.nanstd(data_var_per_film/data_mean_per_film, axis=0),
                                                            'FF UB'      : np.nanmean(data_var_per_film/data_mean_per_film, axis=0) + np.nanstd(data_var_per_film/data_mean_per_film, axis=0),
                                    
                                    
                                                            'COV'        : np.nanmean(np.sqrt(data_var_per_film)/data_mean_per_film, axis=0),
                                                            'COV LB'     : np.nanmean(np.sqrt(data_var_per_film)/data_mean_per_film, axis=0) - np.nanstd(np.sqrt(data_var_per_film)/data_mean_per_film, axis=0),
                                                            'COV UB'     : np.nanmean(np.sqrt(data_var_per_film)/data_mean_per_film, axis=0) + np.nanstd(np.sqrt(data_var_per_film)/data_mean_per_film, axis=0),
                                       
                                                            'CV**2'      : np.nanmean(data_var_per_film/data_mean_per_film**2, axis=0),
                                                            'CV**2 LB'   : np.nanmean(data_var_per_film/data_mean_per_film**2, axis=0) - np.nanstd(data_var_per_film/data_mean_per_film**2, axis=0),
                                                            'CV**2 UB'   : np.nanmean(data_var_per_film/data_mean_per_film**2, axis=0) + np.nanstd(data_var_per_film/data_mean_per_film**2, axis=0),
                                                            
                                                            })
            
            df.to_excel(writer, sheet_name= 'Averaged', index=False, startrow = 0, startcol = 0)
            df.to_excel(writer_all, sheet_name= 'Averaged ' + lgd_name_list[pheno_i], index=False, startrow = 0, startcol = 0)
            """ movie per movie """
            for movie_i in range(len(lgd_per_movie_list[pheno_i])):
                
                df                              = pd.DataFrame({'Time'       : ['T = '+ str(i) + ' min' for i in time_points_list[pheno_i]],
                                                                'Mean'       : np.nanmean(data_per_movie[pheno_i][movie_i], axis=1),
                                                                'Variance'   : np.nanvar(data_per_movie[pheno_i][movie_i], axis=1),
                                                                'FF'         : np.nanvar(data_per_movie[pheno_i][movie_i], axis=1)/np.nanmean(data_per_movie[pheno_i][movie_i], axis=1),
                                                                'COV'        : np.nanstd(data_per_movie[pheno_i][movie_i], axis=1)/np.nanmean(data_per_movie[pheno_i][movie_i], axis=1),
                                                                'CV**2'      : np.nanvar(data_per_movie[pheno_i][movie_i], axis=1)/np.nanmean(data_per_movie[pheno_i][movie_i], axis=1)**2,
                                                                })
        
                df_data                         = pd.DataFrame(data_per_movie[pheno_i][movie_i])
                df_data.columns                 = ['N ' + str(i) for i in range(data_per_movie[pheno_i][movie_i].shape[1])]
                
                df                              = pd.concat([df,df_data], axis=1)
                df.to_excel(writer, sheet_name= lgd_per_movie_list[pheno_i][movie_i][:31], index=False, startrow = 0, startcol = 0)
                
            writer.close()
        writer_all.close()
            
    
##################################################################################





class TranscriptionalSignalInNuclei():

    def __init__(self, data_per_pheno_mitosis, data_per_movie_mitosis, output_folder, inputpath, data_used, lgd_name_list, phenotype_gene, lgd_per_movie_list, FrameLen_list):
        self.output_folder = output_folder
        self.inputpath  = inputpath 
        self.outputPath_extension = data_used + '_by_nuclei/'
        self.data_used = data_used
        
        """ create outputFolder """
        outputPath = self.output_folder + self.outputPath_extension 
        if not os.path.exists(outputPath): 
            os.mkdir(outputPath)



    def Organize_Noise_per_Nuclei(self, data_mitosis, name_list, FrameLen_list):
        
        df_data  =  pd.DataFrame({
                                    'Gene name'            : sum([[name_list[gene_i]]*data_mitosis[gene_i].shape[1] for gene_i in range(len(name_list))], []),
                                    'Mean'                 : np.hstack([np.nanmean(data_mitosis[gene_i], axis =0) for gene_i in range(len(name_list))]),
                                    'Var'                  : np.hstack([np.nanvar(data_mitosis[gene_i], axis = 0 ) for gene_i in range(len(name_list))]),
                                    'FF'                   : np.hstack([np.array(np.nanvar(data_mitosis[gene_i], axis = 0))/np.array(np.nanmean(data_mitosis[gene_i], axis =0) ) for gene_i in range(len(name_list))]),
                                    'CV**2'                : np.hstack([np.array(np.nanvar(data_mitosis[gene_i], axis = 0))/np.array(np.nanmean(data_mitosis[gene_i], axis =0))**2 for gene_i in range(len(name_list))]),
                                    'COV'                  : np.hstack([np.sqrt(np.array(np.nanvar(data_mitosis[gene_i], axis = 0)))/np.array(np.nanmean(data_mitosis[gene_i], axis =0) ) for gene_i in range(len(name_list))]),
                                    'AUC'                  : np.hstack([np.trapz(data_mitosis[gene_i], np.arange(0, data_mitosis[gene_i].shape[0])*FrameLen_list[gene_i], axis=0) for gene_i in range(len(name_list))])})

        return df_data
                        



    def Save_NoisePernuclei(self, df_data_all,  data_combined = '', phenotype_i = '', FrameLen = ''):
        if phenotype_i == '':
            xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + '_Info_per_nuclei.xlsx'
        
        else:
            
            xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + '_' + phenotype_i.replace('/', '') + '_Info_per_nuclei.xlsx'
        
        writer              = pd.ExcelWriter(xlsfilename, engine='xlsxwriter') 
        df_data_all.to_excel(writer, sheet_name= 'All data point') 
        writer.sheets['All data point'].set_column(0,7, len("Mean")+2)
        
        ## add averaged 
        files_name          =  list(set(df_data_all['Gene name'].tolist()))
        mean_per_embryo     = np.zeros((len(files_name), )) 
        var_per_embryo      = np.zeros((len(files_name), )) 
        FF_per_embryo       = np.zeros((len(files_name), )) 
        COV_per_embryo      = np.zeros((len(files_name), )) 
        CV2_per_embryo      = np.zeros((len(files_name), )) 
        
        for i in range(len(files_name)):
            df_data_i           = df_data_all[df_data_all['Gene name'] == files_name[i]]
            mean_per_embryo[i]  = np.nanmean(df_data_i['Mean'].values)
            var_per_embryo[i]   = np.nanmean(df_data_i['Var'].values)
            FF_per_embryo[i]    = np.nanmean(df_data_i['Var'].values/df_data_i['Mean'].values)
            COV_per_embryo[i]   = np.nanmean(np.sqrt(df_data_i['Var'].values)/df_data_i['Mean'].values)
            CV2_per_embryo[i]   = np.nanmean(df_data_i['Var'].values/df_data_i['Mean'].values**2)


        df_data_per_embryo = pd.DataFrame({
                            'Gene name' : files_name,
                            'Mean'      : mean_per_embryo,
                            'Var'       : var_per_embryo, 
                            'FF'        : FF_per_embryo,
                            'CV**2'     : CV2_per_embryo,
                            'COV'       : COV_per_embryo,
                            })
            
        df_data_per_embryo.to_excel(writer, sheet_name= 'Averaged') 
        writer.sheets['Averaged'].set_column(0,7, len("Gene name")+2)
        
        df_data_averaged     = pd.DataFrame({'Gene name' : ['Averaged'],
                                    'Mean': np.nanmean(df_data_per_embryo['Mean']),
                                    'Var'   : np.nanmean(df_data_per_embryo['Var']),
                                    'FF'    : np.nanmean(df_data_per_embryo['FF']),
                                    'CV**2'    : np.nanmean(df_data_per_embryo['CV**2']),
                                    'COV'    : np.nanmean(df_data_per_embryo['COV'])}) 
        
        df_data_per_embryo = pd.concat([df_data_per_embryo, df_data_averaged], ignore_index=True)
        
        if data_combined.size > 0 :
                
            df_data_combined     = pd.DataFrame({'Gene name' : ['Combined'],
                                                'Mean': np.nanmean(np.nanmean(data_combined, axis=0)),
                                                'Var'   : np.nanmean(np.nanvar(data_combined, axis=0)),
                                                'FF'    : np.nanmean(np.nanvar(data_combined, axis=0)/np.nanmean(data_combined, axis=0)),
                                                'CV**2'    : np.nanmean(np.nanvar(data_combined, axis=0)/np.nanmean(data_combined, axis=0)**2),
                                                'COV'    : np.nanmean(np.nanstd(data_combined, axis=0)/np.nanmean(data_combined, axis=0))}) 
                
            df_data_per_embryo = pd.concat([df_data_per_embryo, df_data_combined], ignore_index=True)        
        
        df_data_per_embryo.to_excel(writer, sheet_name= 'Averaged') 
        writer.sheets['Averaged'].set_column(0,7, len("Mean")+2)

        ## add AUC analysis
        files_name          =  list(set(df_data_all['Gene name'].tolist()))
        mean_AUC_embryo     = np.zeros((len(files_name), )) 
        var_AUC_embryo      = np.zeros((len(files_name), )) 
        FF_AUC_embryo       = np.zeros((len(files_name), )) 
        COV_AUC_embryo      = np.zeros((len(files_name), )) 
        CV2_AUC_embryo      = np.zeros((len(files_name), )) 
        
        for i in range(len(files_name)):
            df_data_i           = df_data_all[df_data_all['Gene name'] == files_name[i]]
            mean_AUC_embryo[i]  = np.nanmean(df_data_i['AUC'].values)
            var_AUC_embryo[i]   = np.nanvar(df_data_i['AUC'].values)
            FF_AUC_embryo[i]    = np.nanvar(df_data_i['AUC'].values)/np.nanmean(df_data_i['AUC'].values)
            COV_AUC_embryo[i]   = np.sqrt(np.nanvar(df_data_i['AUC'].values))/np.nanmean(df_data_i['AUC'].values)
            CV2_AUC_embryo[i]   = np.nanvar(df_data_i['AUC'].values)/np.nanmean(df_data_i['AUC'].values)**2


        df_data_AUC_embryo = pd.DataFrame({
                            'Gene name' : files_name,
                            'Mean'      : mean_AUC_embryo,
                            'Var'       : var_AUC_embryo, 
                            'FF'        : FF_AUC_embryo,
                            'CV**2'     : CV2_AUC_embryo,
                            'COV'       : COV_AUC_embryo,
                            })


        df_data_averaged     = pd.DataFrame({'Gene name' : ['Averaged'],
                                            'Mean': np.nanmean(df_data_AUC_embryo['Mean']),
                                            'Var'   : np.nanmean(df_data_AUC_embryo['Var']),
                                            'FF'    : np.nanmean(df_data_AUC_embryo['FF']),
                                            'CV**2'    : np.nanmean(df_data_AUC_embryo['CV**2']),
                                            'COV'    : np.nanmean(df_data_AUC_embryo['COV'])}) 
        if data_combined.size > 0 :
            combined_AUC = np.trapz(data_combined, np.arange(0, data_combined.shape[0])*FrameLen, axis=0)
            df_data_combined     = pd.DataFrame({'Gene name' : ['Combined'],
                                                'Mean': np.nanmean(combined_AUC),
                                                'Var'   : np.nanvar(combined_AUC),
                                                'FF'    : np.nanvar(combined_AUC)/np.nanmean(combined_AUC),
                                                'CV**2'    : np.nanvar(combined_AUC)/np.nanmean(combined_AUC)**2,
                                                'COV'    : np.nanstd(combined_AUC)/np.nanmean(combined_AUC)}) 
                
            df_data_AUC_embryo = pd.concat([df_data_AUC_embryo, df_data_combined], ignore_index=True)   
            
        df_data_AUC_embryo = pd.concat([df_data_AUC_embryo, df_data_averaged], ignore_index=True)
        df_data_AUC_embryo.to_excel(writer, sheet_name= 'AUC Averaged') 
        writer.sheets['AUC Averaged'].set_column(0,7, len("Gene name")+2)
        
        writer.close()
        
        
    
    def noise_datapred_total(self, data_per_movie_mitosis, phenotype_gene, lgd_name_list):

        xlsfilename         = self.output_folder + self.outputPath_extension + '/AUC_' + self.data_used + '_total_noise_Info_per_nuclei.xlsx'
        writer              = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')     
        
        AUC = []
        for i in range(len(lgd_name_list)):
            xlsfilename             = self.output_folder + self.outputPath_extension + '/' + self.data_used + "_" + phenotype_gene[i].replace('/', '') + '_Info_per_nuclei.xlsx'
            df_data_i               = pd.read_excel(xlsfilename, sheet_name = 'AUC Averaged')
            df_data_i               = pd.DataFrame(df_data_i.iloc[:,1:])
            
            df_data_i_AUC               = pd.read_excel(xlsfilename, sheet_name = 'All data point')
            AUC.append(df_data_i_AUC['AUC'].values)
            df_data_i.to_excel(writer, sheet_name= lgd_name_list[i]) 
       
        df_AUC               = pd.DataFrame({'Gene name':sum([[lgd_name_list[gene_i]]*len(AUC[gene_i]) for gene_i in range(len(lgd_name_list))], []), 
                                                    'AUC': np.hstack(AUC)})
        df_AUC.to_excel(writer, sheet_name= 'All data point') 
        writer.close()
        
        
        xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + '_total_noise_Info_per_nuclei.xlsx'
        writer              = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')     
        for i in range(len(lgd_name_list)):
            xlsfilename             = self.output_folder + self.outputPath_extension + '/' + self.data_used + "_" + phenotype_gene[i].replace('/', '') + '_Info_per_nuclei.xlsx'
            df_data_i               = pd.read_excel(xlsfilename, sheet_name = 'Averaged')
            df_data_i               = pd.DataFrame(df_data_i.iloc[:,1:])
            df_data_i.to_excel(writer, sheet_name=  lgd_name_list[i])  
        writer.close()            
        

    
    def Plot_NoisePerNuclei(self, df_data, lgd_pheno = '', phenotype_i = ''):
        sns.set(style = 'whitegrid')
        if 'AUC' in df_data.columns:
                
            """ AUC """ 
            h_AUC, axes         = plt.subplots() 
            sns.violinplot(x ='Gene name', y ='AUC', data = df_data)
            plt.subplots_adjust()
            ylim_min, ylim_max  = axes.axes.get_ylim()
            ylim_min            = ylim_min-10
            sns.swarmplot(x ='Gene name', y ='AUC', data = df_data,color= "black")
            axes.set_ylim((ylim_min, ylim_max))
            plt.xticks(fontsize=7)
            plt.title(lgd_pheno)
            plt.tight_layout()
            if lgd_pheno == '':
                h_AUC.savefig(self.output_folder + self.outputPath_extension + 'AUC_PerNuclei.pdf')
            else:
                h_AUC.savefig(self.output_folder + self.outputPath_extension + 'AUC_' + lgd_pheno + '_PerNuclei.pdf')
        
        
    
        """ FF: Fano Factor """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='FF', data = df_data)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='FF', data = df_data,color= "black")
        ylim_min2, ylim_max2  = axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.title(lgd_pheno)
        plt.tight_layout()
        if lgd_pheno == '':
            h.savefig(self.output_folder + self.outputPath_extension + 'FF_PerNuclei.pdf')
        else:
            h.savefig(self.output_folder + self.outputPath_extension + 'FF_' + lgd_pheno + '_PerNuclei.pdf')
        
        
        """ COV """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='COV', data = df_data)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='COV', data = df_data,color= "black")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.title(lgd_pheno)
        plt.tight_layout()
        if lgd_pheno == '':
            h.savefig(self.output_folder + self.outputPath_extension  + 'COV_PerNuclei.pdf')
        else:
            h.savefig(self.output_folder + self.outputPath_extension + 'COV_' + lgd_pheno + '_PerNuclei.pdf')
        
        
        
        """ CV**2 """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='CV**2', data = df_data)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='CV**2', data = df_data,color= "black")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.title(lgd_pheno)
        plt.tight_layout()
        if lgd_pheno == '':
            h.savefig(self.output_folder + self.outputPath_extension + 'CV_PerNuclei.pdf')
        else:
            h.savefig(self.output_folder + self.outputPath_extension + 'CV_' + lgd_pheno + '_PerNuclei.pdf')
        
        
        
        
        """ mean Factor """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='Mean', data = df_data)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='Mean', data = df_data,color= "black")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.title(lgd_pheno)
        plt.tight_layout()
        if lgd_pheno == '':
            h.savefig(self.output_folder + self.outputPath_extension + 'Mean_PerNuclei.pdf')
        else:
            h.savefig(self.output_folder + self.outputPath_extension + 'Mean_' + lgd_pheno + '_PerNuclei.pdf')
        
        
        

    

        
        """  VAriance """ 
        h, axes             = plt.subplots() 
        sns.violinplot(x ='Gene name', y ='Var', data = df_data)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='Var', data = df_data,color= "black")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.title(lgd_pheno)
        plt.tight_layout()
        if lgd_pheno == '':
            h.savefig(self.output_folder + self.outputPath_extension + 'Var_PerNuclei.pdf')
        else:
            h.savefig(self.output_folder + self.outputPath_extension + 'Var_' + lgd_pheno + '_PerNuclei.pdf')
        
        
        
        
    def plotting_results_total_violonPlot(self, filename, extension = ''):
        xls = pd.ExcelFile(filename)
        
        """ plotting all Pol II """
        if 'All data point' in xls.sheet_names:
            df_data = pd.read_excel(filename, sheet_name = 'All data point')
            if 'Nbr Pol_II' in df_data.columns:
                
                sns.set(style = 'whitegrid')
        
                h, axes         = plt.subplots() 
                sns.violinplot(x ='Gene name', y ='Nbr Pol_II', data = df_data)
                plt.subplots_adjust()
                ylim_min, ylim_max  = axes.axes.get_ylim()
                ylim_min            = ylim_min-2
                sns.swarmplot(x ='Gene name', y ='Nbr Pol_II', data = df_data,color= "black")
                axes.set_ylim((ylim_min, ylim_max))
                plt.xticks(fontsize=7)
                plt.tight_layout()
                h.savefig(self.output_folder + self.outputPath_extension + 'polII_PerNuclei.pdf')
        
            else:
            
                sns.set(style = 'whitegrid')
             
                h, axes         = plt.subplots() 
                sns.violinplot(x ='Gene name', y ='AUC', data = df_data)
                plt.subplots_adjust()
                ylim_min, ylim_max  = axes.axes.get_ylim()
                ylim_min            = ylim_min-2
                sns.swarmplot(x ='Gene name', y ='AUC', data = df_data,color= "black")
                axes.set_ylim((ylim_min, ylim_max))
                plt.xticks(fontsize=7)
                plt.tight_layout()
                h.savefig(self.output_folder + self.outputPath_extension + 'AUC_PerNuclei.pdf')           
        
        sheet_names = xls.sheet_names
        sheet_names = [name for name in sheet_names if name != 'All data point']
        
        name_list = []
        
        df_Average = pd.DataFrame({'Mean':[],	'Var':[],	'FF':[],	'COV':[],	'CV**2':[]})
        df_Combined = pd.DataFrame({'Mean':[],	'Var':[],	'FF':[],	'COV':[],	'CV**2':[]})
        df = pd.DataFrame({'Mean':[],	'Var':[],	'FF':[],	'COV':[],	'CV**2':[]})
        
        for i in range(len(sheet_names)):
            df_i = pd.read_excel(filename, sheet_name= sheet_names[i])
            df = pd.concat([df, df_i.iloc[:-2,:]], ignore_index=True)
            
            df_Average = pd.concat([df_Average, df_i.iloc[-2,:].to_frame().T], ignore_index=True)
            df_Combined = pd.concat([df_Combined, df_i.iloc[-1,:].to_frame().T], ignore_index=True)
            
            name_list = name_list+ [sheet_names[i]]*df_i.iloc[:-2,:].shape[0]
        
        df_Average['Gene name'] = sheet_names
        df_Combined['Gene name'] = sheet_names
        df['Gene name'] = name_list
    
        sns.set(style = 'whitegrid')
        
        
        """ FF: Fano Factor """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='FF', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='FF', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='FF', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='FF', data = df_Combined,color= "red")

        ylim_min2, ylim_max2  = axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'FF_PerMovie.pdf')

        
        """ COV """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='COV', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='COV', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='COV', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='COV', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'COV_PerMovie.pdf')
        
        
        
        """ CV**2 """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='CV**2', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='CV**2', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='CV**2', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='CV**2', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'CV_PerMovie.pdf')
        
        
        
        """ mean Factor """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='Mean', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='Mean', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='Mean', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='Mean', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'Mean_PerMovie.pdf')
        
        
        
        """  VAriance """ 
        h, axes             = plt.subplots() 
        sns.violinplot(x ='Gene name', y ='Var', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='Var', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='Var', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='Var', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'Var_PerMovie.pdf')
        
        
            
        
class PosPredInNuclei():
    def __init__(self, data_per_pheno_mitosis, data_per_movie_mitosis, output_folder, inputpath, data_used, lgd_name_list, phenotype_gene, lgd_per_movie_list, FrameLen_list):

        self.output_folder          =  output_folder
        self.inputpath              =  inputpath
        self.outputPath_extension   =  data_used + '_by_nuclei/'
        self.data_used              =  data_used
        
        """ create outputFolder """
        outputPath = self.output_folder + self.outputPath_extension 
        if not os.path.exists(outputPath): 
            os.mkdir(outputPath)


    
    def Organize_Noise_per_Embryo_Pol_II(self, data_mitosis, data_combined, name_list, FrameLen_list, phenotype_i = ''):
        if phenotype_i == '':
            xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + '_Info_per_nuclei.xlsx'
        
        else:
            
            xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + "_" + phenotype_i.replace('/', '') + '_Info_per_nuclei.xlsx'

        writer              = pd.ExcelWriter(xlsfilename, engine='xlsxwriter') 
        
        for gene_i in range(len(name_list)):
            df_data_per_embryo = pd.DataFrame({
                                'Gene name'  : [' N ' + str(i) for i in range(data_mitosis[gene_i].shape[0])],
                                'Nbr Pol II' : data_mitosis[gene_i],
                                })
            df_data_per_embryo.to_excel(writer, sheet_name= name_list[gene_i][:31]) 
            
        df_data             = pd.DataFrame({
                            'Gene name'         : name_list ,
                            'Mean'              : np.hstack([np.nanmean(data_mitosis[gene_i], axis=0) for gene_i in range(len(name_list))]),
                            'Var'               : np.hstack([np.nanvar(data_mitosis[gene_i], axis=0) for gene_i in range(len(name_list))]),
                            'FF'                : np.hstack([np.array(np.nanvar(data_mitosis[gene_i], axis=0)) / np.array(np.nanmean(data_mitosis[gene_i], axis=0)) for gene_i in range(len(name_list))]),
                            'CV**2'             : np.hstack([np.array(np.nanvar(data_mitosis[gene_i], axis=0)) / np.array(np.nanmean(data_mitosis[gene_i], axis=0))**2 for gene_i in range(len(name_list))]),
                            'COV'               : np.hstack([np.sqrt(np.array(np.nanvar(data_mitosis[gene_i], axis=0))) / np.array(np.nanmean(data_mitosis[gene_i], axis=0)) for gene_i in range(len(name_list))])
                            })

        df_data_averaged    = pd.DataFrame({'Gene name'  : ['Averaged'],
                                            'Mean'       : np.nanmean(df_data['Mean']),
                                            'Var'        : np.nanmean(df_data['Var']),
                                            'FF'         : np.nanmean(df_data['FF']),
                                            'CV**2'      : np.nanmean(df_data['CV**2']),
                                            'COV'        : np.nanmean(df_data['COV'])})

        df_data = pd.concat([df_data, df_data_averaged], ignore_index=True)

        df_data_combined     = pd.DataFrame({'Gene name' : ['Combined'],
                                            'Mean'  : np.nanmean(data_combined),
                                            'Var'   : np.nanvar(data_combined),
                                            'FF'    : np.nanvar(data_combined) / np.nanmean(data_combined),
                                            'CV**2' : np.nanvar(data_combined) / np.nanmean(data_combined)**2,
                                            'COV'   : np.nanstd(data_combined) / np.nanmean(data_combined)})

        df_data  =  pd.concat([df_data, df_data_combined], ignore_index=True)

        df_data.to_excel(writer, sheet_name='Averaged')
        writer.sheets['Averaged'].set_column(0, 7, len("Mean") + 2)

        df_data = pd.DataFrame({'Gene name': sum([[name_list[gene_i] for i in range(data_mitosis[gene_i].shape[0])] for gene_i in range(len(name_list))], []), 'Nbr Pol_II': np.hstack(data_mitosis)})
        df_data.to_excel(writer, sheet_name= 'All data point') 
        writer.sheets['All data point'].set_column(0,7, len("Mean") + 2)
        writer.close()

        return df_data

    def noise_pospres_total(self, pol_II_nbr_TW_per_movie, phenotype_gene, lgd_name_list):
        xlsfilename  =  self.output_folder + self.outputPath_extension + '/' + self.data_used + '_total_noise_Info_per_nuclei.xlsx'
        writer       =  pd.ExcelWriter(xlsfilename, engine='xlsxwriter')

        df_results   =  pd.DataFrame({ 'Gene name': sum([[lgd_name_list[gene_i]]*len(np.hstack(pol_II_nbr_TW_per_movie[gene_i])) for gene_i in range(len(lgd_name_list))], []),
                        'Nbr Pol_II': sum([np.hstack(pol_II_nbr_TW_per_movie[gene_i]).tolist() for gene_i in range(len(lgd_name_list))], [])})
        df_results.to_excel(writer, sheet_name= 'All data point') 

        for i in range(len(lgd_name_list)):
            xlsfilename             = self.output_folder + self.outputPath_extension + '/' + self.data_used + "_" +phenotype_gene[i].replace('/', '') + '_Info_per_nuclei.xlsx'
            df_data_i               = pd.read_excel(xlsfilename, sheet_name = 'Averaged')
            df_data_i               = pd.DataFrame(df_data_i.iloc[:,1:])
            df_data_i.to_excel(writer, sheet_name=lgd_name_list[i])

        writer.close()
        
        
        
    def Plot_NoisePerNuclei_Pol_II(self, df_data, lgd_pheno = '', phenotype_i = ''):
        if phenotype_i == '':
            xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + '_Info_per_nuclei.xlsx'

        else:
            
            xlsfilename         = self.output_folder + self.outputPath_extension + '/' + self.data_used + "_"+  phenotype_i.replace('/', '') + '_Info_per_nuclei.xlsx'


        df_data = pd.read_excel(xlsfilename, sheet_name = 'All data point')
        sns.set(style = 'whitegrid')

        h, axes         = plt.subplots() 
        sns.violinplot(x ='Gene name', y ='Nbr Pol_II', data = df_data)
        plt.subplots_adjust()
        ylim_min, ylim_max  = axes.axes.get_ylim()
        ylim_min            = ylim_min-2
        sns.swarmplot(x ='Gene name', y ='Nbr Pol_II', data = df_data,color= "black")
        axes.set_ylim((ylim_min, ylim_max))
        plt.xticks(fontsize=7)
        plt.title(lgd_pheno)
        plt.tight_layout()
        if lgd_pheno == '':
            h.savefig(self.output_folder + self.outputPath_extension + 'polII_PerNuclei.pdf')
        else:
            h.savefig(self.output_folder + self.outputPath_extension +  lgd_pheno + '_polII_PerNuclei.pdf')
            
            
            
            
        df_data = pd.read_excel(xlsfilename, sheet_name = 'Averaged')
        
        # Plot the grouped bar chart
        fig, ax = plt.subplots(2,3, figsize=(10, 4))
        df_results_mean = df_data[['Mean']].T
        df_results_mean.plot(kind='bar', rot=0, width=0.7, ax=ax[0,0])
        ax[0,0]
        ax[0,0].set_ylabel('Values')
        ax[0,0].legend_.remove()
        
        df_results_var = df_data[['Var']].T
        df_results_var.plot(kind='bar', rot=0, width=0.7, ax = ax[0,1])
        ax[0,1].legend_.remove()
        
        df_results_FF = df_data[['FF']].T
        df_results_FF.plot(kind='bar', rot=0, width=0.7, ax = ax[1,0])
        ax[1,0].set_ylabel('Values')
        ax[1,0].legend_.remove()
        
        df_results_COV = df_data[['COV']].T
        df_results_COV.plot(kind='bar', rot=0, width=0.7, ax = ax[1,1])
        ax[1,1].legend_.remove()
        
        df_results_CV2 = df_data[['CV**2']].T
        bar_container = df_results_CV2.plot(kind='bar', rot=0, width=0.7, ax= ax[1,2])
        ax[1,2].legend_.remove()
        
        color_container = bar_container.containers
        colors = []
        for i in range(len(color_container)):
            for bar in bar_container.containers[i]: 
                colors.append(bar.get_facecolor())
                
        # Define custom legend labels and colors
        legend_elements = []
        for i in range(len(color_container)):
            legend_elements.append(Patch(facecolor=colors[i], label=list(df_results_CV2.columns)[i]))
        
        ax[0, 2].legend(handles=legend_elements, loc='center', prop={'size': 8})
        ax[0, 2].set_xticks([])
        ax[0, 2].set_yticks([])
        ax[0, 2].spines['top'].set_visible(False)
        ax[0, 2].spines['right'].set_visible(False)
        ax[0, 2].spines['bottom'].set_visible(False)
        ax[0, 2].spines['left'].set_visible(False)
        ax[0, 2].grid(False)
        
        
        plt.tight_layout()
        fig.savefig(self.output_folder + self.outputPath_extension + 'PolII_all_results_' + lgd_pheno + '.pdf')



        
    
    def plotting_results_total_violonPlot(self, filename, extension = ''):
        xls = pd.ExcelFile(filename)
        
        """ plotting all Pol II """
        if 'All data point' in xls.sheet_names:
            
            df_data = pd.read_excel(filename, sheet_name = 'All data point')
            sns.set(style = 'whitegrid')
    
            h, axes         = plt.subplots() 
            sns.violinplot(x ='Gene name', y ='Nbr Pol_II', data = df_data)
            plt.subplots_adjust()
            ylim_min, ylim_max  = axes.axes.get_ylim()
            ylim_min            = ylim_min-2
            sns.swarmplot(x ='Gene name', y ='Nbr Pol_II', data = df_data,color= "black")
            axes.set_ylim((ylim_min, ylim_max))
            plt.xticks(fontsize=7)
            plt.tight_layout()
            h.savefig(self.output_folder + self.outputPath_extension + 'polII_PerNuclei.pdf')
    
            
        
        sheet_names = xls.sheet_names
        sheet_names = [name for name in sheet_names if name != 'All data point']
        
        name_list = []
        
        df_Average = pd.DataFrame({'Mean':[],	'Var':[],	'FF':[],	'COV':[],	'CV**2':[]})
        df_Combined = pd.DataFrame({'Mean':[],	'Var':[],	'FF':[],	'COV':[],	'CV**2':[]})
        df = pd.DataFrame({'Mean':[],	'Var':[],	'FF':[],	'COV':[],	'CV**2':[]})
        
        for i in range(len(sheet_names)):
            df_i = pd.read_excel(filename, sheet_name= sheet_names[i])
            df = pd.concat([df, df_i.iloc[:-2,:]], ignore_index=True)
            
            df_Average = pd.concat([df_Average, df_i.iloc[-2,:].to_frame().T], ignore_index=True)
            df_Combined = pd.concat([df_Combined, df_i.iloc[-1,:].to_frame().T], ignore_index=True)
            
            name_list = name_list+ [sheet_names[i]]*df_i.iloc[:-2,:].shape[0]
        
        df_Average['Gene name'] = sheet_names
        df_Combined['Gene name'] = sheet_names
        df['Gene name'] = name_list
    
        sns.set(style = 'whitegrid')
        
        
        """ FF: Fano Factor """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='FF', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='FF', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='FF', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='FF', data = df_Combined,color= "red")

        ylim_min2, ylim_max2  = axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'FF_PerMovie.pdf')

        
        """ COV """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='COV', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='COV', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='COV', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='COV', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'COV_PerMovie.pdf')
        
        
        
        """ CV**2 """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='CV**2', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='CV**2', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='CV**2', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='CV**2', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'CV_PerMovie.pdf')
        
        
        
        """ mean Factor """ 
        h, axes             = plt.subplots()
        
        sns.violinplot(x ='Gene name', y ='Mean', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='Mean', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='Mean', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='Mean', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'Mean_PerMovie.pdf')
        
        
        
        """  VAriance """ 
        h, axes             = plt.subplots() 
        sns.violinplot(x ='Gene name', y ='Var', data = df)
        plt.subplots_adjust()
        ylim_min1, ylim_max1  = axes.axes.get_ylim()

        sns.swarmplot(x ='Gene name', y ='Var', data = df,color= "black")
        sns.swarmplot(x ='Gene name', y ='Var', data = df_Average,color= "white")
        sns.swarmplot(x ='Gene name', y ='Var', data = df_Combined,color= "red")
        ylim_min2, ylim_max2  = axes.axes.get_ylim()
        
        axes.set_ylim((min(ylim_min1, ylim_min2), max(ylim_max1, ylim_max2)))
        plt.xticks(fontsize=7)
        plt.tight_layout()
        h.savefig(self.output_folder + self.outputPath_extension + extension + 'Var_PerMovie.pdf')
        
        
    



        
    


    
    
        
########## 
        
class PosPredPreparation():
    def __init__(self, inputpath, extension, data_type, phenotype_gene, phenotype_rawData, phenotype_BurstDeconv, fParam_name, tfinal, TW_min = 1):
        self.TW_min                 = TW_min
        self.inputpath              = inputpath
        self.extension              = extension
        self.data_type              = data_type
                
        

        
        
        data_per_pheno_mitosis = [] # contains list of data averaged in  a time window of TW_min for combined movies w.r.t. mitosis of each phenotype 
        data_per_movie_mitosis = [] # contains list of data averaged in  a time window of TW_min for combined movies w.r.t. 20 % of first activation of each phenotype 
                
        lgd_name_list          = [] #c contains list of the name of the phenotype
        time_points_list       = [] # contains the time points analysed of each phenotype
        lgd_per_movie_list     = [] # contains list of list of the movie name of each phenotype
        FrameLen_list          = []
        

        
        """ movies combining """
        
        name_split             = phenotype_gene[len(phenotype_gene)-1].split('/')
        if len(phenotype_gene) == 1:
            index_legend = 0
        else:
            
            for i in range(len(name_split)-1):
                if phenotype_gene[-1].split('/')[i] != phenotype_gene[-2].split('/')[i]:
                    index_legend = i

            
        for gene_i in range(len(phenotype_BurstDeconv)): 
            fParam                  = self.inputpath + phenotype_BurstDeconv[gene_i] + fParam_name
            parameter_content       = np.load(fParam)
            
            Polym_speed        = parameter_content['Polym_speed']
            DureeSignal        = parameter_content['DureeSignal']
            TaillePostMarq     = parameter_content['TaillePostMarq']
            TailleSeqMarq      = parameter_content['TailleSeqMarq']
            EspaceInterPolyMin = parameter_content['EspaceInterPolyMin']
            FreqEchSimu        = 1/(EspaceInterPolyMin/Polym_speed)
            FreqEchSimu             = parameter_content['FreqEchSimu']
            FrameLen                = parameter_content['FrameLen']
            
            signal_length      = int((TailleSeqMarq+TaillePostMarq)/Polym_speed)
            Tfinal_bp          = int((tfinal*60 + DureeSignal)*FreqEchSimu)

            xlsPath                                                  = self.inputpath + phenotype_rawData[gene_i] 
            respath                                                  = self.inputpath + phenotype_BurstDeconv[gene_i] + 'resultDec' + self.data_type + '/'
            
            filesList                                                = np.array(os.listdir(xlsPath )) # list of the data
            filesList                                                = [filesList[i].replace(self.extension,'') for i in range(len(filesList) )]
        
            [DataExp_combined, DataPred_combined, PosPred_combined, PosPred_combined_corrected, T0_combined, tmax_combined, FrameLen] = movies_combining(xlsPath,respath, filesList, fParam, self.extension)

            data_per_pheno_mitosis.append(PosPred_combined[:Tfinal_bp,:])
            
            time_points                                             = np.arange(0,DataExp_combined.shape[0]*FrameLen,self.TW_min*60)
            FrameLen_list.append(FrameLen)
            
            lgd_name_list.append(phenotype_gene[gene_i].split('/')[index_legend])
            lgd_per_movie_list.append([clean_movie_name(i) for i in filesList])    
        
        
            pol_2_nbr_per_movie_mitosis                             = []
            
            for data_i in range(len(filesList)):       
                resultPath_i                                                            = respath +  'result_' + filesList[data_i].replace('_','') + '.npz'
                xlsFile_i                                                               = xlsPath + filesList[data_i] + self.extension
                [T0, tmax, DataExp_i, DataPred_i, PosPred_i]                            = data_reshape(xlsFile_i, resultPath_i, fParam, 60)
                
                """ movie i, t0 = begining of movie """               
                pol_2_nbr_per_movie_mitosis.append(PosPred_i[:Tfinal_bp,:])
            
                
            data_per_movie_mitosis.append(pol_2_nbr_per_movie_mitosis)
            time_points_list.append((time_points + signal_length)/60)
        
        
        self.data_per_pheno_mitosis = data_per_pheno_mitosis
        self.data_per_movie_mitosis = data_per_movie_mitosis
        
        self.time_points_list = time_points_list
        self.FrameLen_list = FrameLen_list
        
        self.lgd_name_list = lgd_name_list
        self.lgd_per_movie_list = lgd_per_movie_list

        
###############################################################################
        
class TranscriptionSignalPreparation():
    def __init__(self,inputpath, data_used, tfinal, TW_min, phenotype_rawData, phenotype_gene, data_type, extension, calibration=1, color_heatmap = 'YlOrBr', phenotype_BurstDeconv = '', fParam_name = ''):    
        """ define instance variable to store output """
     
        self.data_per_pheno_mitosis = [] # contains list of data averaged in  a time window of TW_min for combined movies w.r.t. mitosis of each phenotype 
        self.data_per_movie_mitosis = [] # contains list of data averaged in  a time window of TW_min for combined movies w.r.t. 20 % of first activation of each phenotype 
        
        
        self.data_per_pheno_aligned = []  # contains list of list of data averaged in  a time window of TW_min for each movies w.r.t. mitosis for each phenotype: data[i] contains a list of each movie of phenotype i
        self.data_per_movie_aligned = []  # contains list of list of data averaged in  a time window of TW_min for each movies w.r.t. 20% of first activation for each phenotype: data[i] contains a list of each movie of phenotype i
        
        
        self.lgd_name_list          = [] #c contains list of the name of the phenotype
        self.FrameLen_list          = [] # contains the framelen of each phenotype
        self.lgd_per_movie_list     = [] # contains list of list of the movie name of each phenotype
        self.time_points_list       = [] # contains the time points of each phenotype
        
        name_split = phenotype_gene[len(phenotype_gene)-1].split('/')
        if len(phenotype_gene) == 1:
            index_legend = 0
        else:
            
            for i in range(len(name_split)-1):
                if phenotype_gene[-2].split('/')[i] != phenotype_gene[-1].split('/')[i]:
                    index_legend = i
            
            
        for gene_i in range(len(phenotype_gene)):
            
            lgd_name                                                 = phenotype_gene[gene_i].split('/')[index_legend]
            
            xlsPath                                                  = inputpath + phenotype_rawData[gene_i] 
                
            filesList                                                = np.array(os.listdir(xlsPath )) # list of the data
            filesList                                                = [filesList[i].replace(extension,'') for i in range(len(filesList) )]
            self.lgd_per_movie_list.append([clean_movie_name(i) for i in filesList])    
            
            """ combine movies"""
            if data_used == 'DataExp':
                [FrameLen, data_mitosis_combined, tmax_combined, tstart]                                                                      = movies_combining_rawData(xlsPath, calibration, extension)
            else:
                fParam                                                                                                                        = inputpath + phenotype_BurstDeconv[gene_i] + fParam_name
                respath                                                                                                                       = inputpath + phenotype_BurstDeconv[gene_i] + 'resultDec' + data_type + '/'
                [DataExp_combined, data_mitosis_combined, PosPred_combined, PosPred_combined_corrected, T0_combined, tmax_combined, FrameLen] = movies_combining(xlsPath,respath, filesList, fParam, extension)
                
                parameter_content       = np.load(fParam)
            
                Polym_speed        = parameter_content['Polym_speed']
                TaillePostMarq     = parameter_content['TaillePostMarq']
                TailleSeqMarq      = parameter_content['TailleSeqMarq']
                signal_length      = int((TailleSeqMarq+TaillePostMarq)/Polym_speed)     
                time_points                                             = np.arange(0,DataExp_combined.shape[0]*FrameLen, TW_min*60)
            
            """ compute the data in a TW for combined movies w.r.t. mitosis """
            Tfinal_bp                                               = int(tfinal*60/FrameLen)
            
            
            self.data_per_pheno_mitosis.append(data_mitosis_combined[:Tfinal_bp,:]) 
            self.lgd_name_list.append(lgd_name)
            self.FrameLen_list.append(FrameLen)
            if data_used != 'DataExp':
                self.time_points_list.append((time_points + signal_length)/60)
            
            """ compute the data in a time windows for a combined data w.r.t. 1st activation """
            data_aligned_combined                                    = allign_data_Activation(data_mitosis_combined[:Tfinal_bp,:]) # allign the data
            self.data_per_pheno_aligned.append(data_aligned_combined)      
            
            
            """ compute the data in a time windows for mitosis and w.r.t. 1 st activation movie per movie """
            
            data_i_TW_perPheno_mitosis                               = []
            data_i_TW_perPheno_aligned                               = []
            for i in range(len(filesList)):
                xlsFile_i                                                               = xlsPath + filesList[i] + extension                
                # allign w.r.t. mitosis
                if data_used == 'DataExp':
                    [T0_i, tmax_i, data_mitosis_i]                                          = data_reshape_raw(xlsFile_i)
                else:
                    resultPath_i                                                            = respath +  'result_' + filesList[i].replace('_','') + '.npz'
                    [T0, tmax, DataExp, data_mitosis_i, PosPred]                            = data_reshape(xlsFile_i, resultPath_i, fParam)
    
                data_i_TW_perPheno_mitosis.append(data_mitosis_i[:Tfinal_bp,:])
                
                # allign w.r.t. 20% of first activation
                data_aligned_i                                                          = allign_data_Activation(data_mitosis_i[:Tfinal_bp,:]) # allign the data
                data_i_TW_perPheno_aligned.append(data_aligned_i) 
            
            
            self.data_per_movie_mitosis.append(data_i_TW_perPheno_mitosis)
            self.data_per_movie_aligned.append(data_i_TW_perPheno_aligned)
