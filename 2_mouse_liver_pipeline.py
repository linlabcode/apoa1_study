#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2019 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run scripts for processing liver mouse data for Lagor project




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================

py27_path = '/storage/cylin/anaconda3/envs/py27_anaconda/bin/python2'

projectName = 'Lagor/For_Charles/'
genome ='mm10'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files
maskFile = '/storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Annotation/Masks/mm10_blacklist.bed'

#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Seq
chip_data_file = '%sfinal_data_tables/mm10_liver_CHIP_data_table.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)


    print('\n\n')
    print('#======================================================================')
    print('#==========================II. RUNNING MACS============================')
    print('#======================================================================')
    print('\n\n')

    data_dict=  pipeline_dfci.loadDataTable(chip_data_file)

    k27ac_list= [name for name in data_dict.keys() if name.count('27ac') == 1 and name.upper().count('WCE') == 0]
    
    for name in k27ac_list:
        print(data_dict[name])

    
    pipeline_dfci.run_macs(chip_data_file,projectFolder,macsFolder,macsEnrichedFolder,wiggleFolder,True,k27ac_list)


    print('\n\n')
    print('#======================================================================')
    print('#===================III. CALL ROSE INDIVIDUALLY========================')
    print('#======================================================================')
    print('\n\n')
    analysis_name = 'MOUSE_LIVER_H3K27AC'
    parentFolder = utils.formatFolder('%s%s' % (roseFolder,analysis_name),True)
    #pipeline_dfci.callRose2(chip_data_file,macsEnrichedFolder,parentFolder,namesList=k27ac_list,extraMap = [],inputFile='',tss=2500,stitch='',bashFileName ='',mask=maskFile,useBackground=True,py27_path =py27_path)

    #run rose2 wrapper for both
    enhancer_bashFileName,enhancer_region_map_path,names_list = define_enhancer_landscape(projectFolder,pipeline_dir,chip_data_file,analysis_name ,k27ac_list)
    print(enhancer_bashFileName,enhancer_region_map_path,names_list)

    #runs only if no output detected                                                                      
    if not utils.checkOutput(enhancer_region_map_path,0,0):                                               
        print(enhancer_bashFileName)                                                                      
        os.system('bash %s' % (enhancer_bashFileName))   



    sys.exit()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~DEFINING H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(projectFolder,pipeline_dir,data_file,analysis_name = '',names_list = [],stitch='',use_background = True):

    '''
    defines the NB enhancer baseline using H3K27ac chips 
    enhancers defined using auto optimized stitching of nearby regions
    w/ a 2.5kb tss exclusion
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For H3K27AC
    #with TSS exclusion and auto stitching

    dataDict = pipeline_dfci.loadDataTable(data_file)

    if  analysis_name == '':
        #get the name from the data_file
        analysis_name = data_file.split('/')[-1].split('.')[0]
    
    print("RUNNING ANALYSIS: '%s'" % (analysis_name))


    if names_list == []:
        names_list = [name for name in dataDict.keys() if name.upper().count('H3K27AC') == 1]
    print('FOR H3K27AC DATASETS USING:')
    print(names_list)


    bamFileList = [dataDict[name]['bam'] for name in names_list]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in names_list]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in names_list]
    bedString = string.join(bedFileList,',')

    roseFolder = '%smeta_rose/' % (projectFolder)
    roseFolder = utils.formatFolder(roseFolder,True)

    outputFolder = '%s%s/' % (roseFolder,analysis_name)
    bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysis_name)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = '%s %sROSE2_META.py -g %s -i %s -r %s -c %s -o %s -n %s -t 2500 --mask %s' % (py27_path,pipeline_dir,genome,bedString,bamString,controlBamString,outputFolder,analysis_name,maskFile)
    if stitch != '':
        metaRoseCmd += ' -s %s' % (stitch)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()


    #getting the region_map_path as a way to know if it's done
    region_map_path = '%s%s/%s_AllEnhancers.table.txt' % (roseFolder,analysis_name,analysis_name)
    return bashFileName,region_map_path,names_list




#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
