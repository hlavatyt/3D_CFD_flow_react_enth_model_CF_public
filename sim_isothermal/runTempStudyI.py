#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  runTempStudy.py
#  
#  Copyright 2019 Tomáš Hlavatý & Martin Isoz<hlavatyo@vscht.cz>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
####################################
# Run script for isothermal simulations of flow and DCR in macro-scale CF  

import os
import time
import sys

speciesLst  = ['CO','NO','O2','N2O']


outFile     = './outletData.dat'						# outlet file for results
isOnWallIn0org =0 										# 1 - isOnWallCoating and isInWallCoating fields are prepared in the 0.org folder, 0 - not prepared
nIterFlow    = 1500										# number of iterations for the first run of flow simulation
# ~nIterFlow    = 5
nIterFlowThen = 250									# number of iterations for other runs of flow simulations
# ~nIterFlowThen = 4
nIterKin     = 50												#number of iterations for DCR model
# ~nIterKin     = 6

patchName = 'average'                                 # different patch names for postProcesing - for kraken
# ~patchName = 'areaAverage'                    # different patch names for postProcesing - my PC                     

#tempLst with solved temperatures
tempLst = [303.0, 333.0, 363.0, 393.0, 423.0, 453.0, 483.0, 513.0, 543.0, 573.0, 603.0, 633.0, 663.0, 693.0, 723.0, 753.0, 783.0, 813.0, 843.0, 873.0] 	
# ~ tempLst = [513.0]
startTemp = tempLst[0]

runParallel = True											# paralel run, Note: I am using parallel run, not checke if the non-paralallel run works now, earlier worked
#Note: only porousSimpleFoams runs in parallel, DCR not
nProc = 8														
# ~nProc = 16

# auxiliary func for parse of folder names - very ugly but working, see code below for description
def isInt(value):
  try:
    int(value)
    return True
  except ValueError:
    return False    

# -- time the code #####################################################
startTime = time.time()

# -- rm outfile#########################################################
os.system('rm %s'%outFile)

# -- foreach temp#######################################################
for ind2 in range(len(tempLst)):
    temp = tempLst[ind2]
    # kinematic viscosity###############################################
    #nu = (9.57353896655317e-5 * (temp-273.15)**2 + 0.0847135666774343 * (temp-273.15) +13.7018126018104) *10**(-6)
    mu = ((1.458e-6)*(temp)**(1.0/2.0))/(1.0+110.4/temp)
    print mu
    # density of air####################################################
    rho = 101325.0/(287.058*temp)
    nu = mu/rho
    # -- read and update temperature and nu in the transportProperties##
    print ('Updating T in transportProperties to %5.2f and nu to %5.8f'%(temp,nu))
    fileHn  = open('constant/transportProperties','r')
    lineLst = fileHn.readlines()
    for ind in range(len(lineLst)):
        if lineLst[ind].find('cTemp') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            #~ lineLst[ind] =  ' '.join(lineLst[ind].split(' ')[1:len(lineLst[ind].split(' '))]) + ' %5.2f;\n'%temp
            lineLst[ind] = '\tcTemp\tcTemp\t[0 0 0 1 0 0 0]\t%5.2f;\n'%temp
        if lineLst[ind].find('nu\t') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed nu
            lineLst[ind] = 'nu\t[0 2 -1 0 0 0 0]\t%5.8f;\n'%nu
    fileHn.close()
    with open('constant/transportProperties','w') as fileHn:
        fileHn.writelines(lineLst)   
    fileHn.close()
    
    # if first run######################################################
    if startTemp == temp:
        # --clean of old cases##########################################
        os.system('./Allclean')
        
        # -- parallel run ##############################################
        if runParallel:
            print('Running parallel on %d proc'%nProc)
            # --set number of proc######################################
            print ('Updating decomposeParDict numberOfSubdomains to %d'%nProc)
            fileHn  = open('system/decomposeParDict','r')
            lineLst = fileHn.readlines()
            for ind in range(len(lineLst)):
                if lineLst[ind].find('numberOfSubdomains') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
                    lineLst[ind] = 'numberOfSubdomains\t%d;\n'%nProc
            fileHn.close()
            with open('system/decomposeParDict','w') as fileHn:
                fileHn.writelines(lineLst)
            fileHn.close()
            
        # ~# -- checkWriteInterval
        # ~print ('Checking writeInterval in controlDict')
        # ~fileHn  = open('system/controlDict','r')
        # ~lineLst = fileHn.readlines()
        # ~for ind in range(len(lineLst)):
            # ~if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                # ~spLine = lineLst[ind].split('\t')
                # ~intSpLine = int(spLine[1].replace(';\n',''))
                # ~if min(nIterFlow, nIterKin, nIterFlowThen) > float(spline[1]):
                # ~lineLst[ind] = 'writeInterval \t %d;\n' %min(nIterFlow, nIterKin, nIterFlowThen, intSpLine)
                # ~print ('WriteInterval changed to %d' %min(nIterFlow, nIterKin, nIterFlowThen, intSpLine))
        # ~fileHn.close()
        # ~with open('system/controlDict','w') as fileHn:
            # ~fileHn.writelines(lineLst)
        # ~fileHn.close()
        
        # -- set number of iterations of flow in controlDict############
        # ~print ('Updating controlDict endtime to %d'%nIterFlow)
        # ~fileHn  = open('system/controlDict','r')
        # ~lineLst = fileHn.readlines()
        # ~for ind in range(len(lineLst)):
            # ~if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                # ~lineLst[ind] = 'endTime\t%d;\n'%nIterFlow
            # ~if lineLst[ind].find('application') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
                # ~lineLst[ind] = 'application\tporousSimpleFoam;\n'
            # ~if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                # ~spLine = lineLst[ind].split('\t')
                # ~intSpLine = int(spLine[1].replace(';\n',''))
                # ~lineLst[ind] = 'writeInterval \t %d;\n' %nIterFlow
        # ~fileHn.close()
        # ~with open('system/controlDict','w') as fileHn:
            # ~fileHn.writelines(lineLst)
        # ~fileHn.close()
        
        # -- running blockMesh, stitchMesh, createPatch etc.############
        print ('Running blockMesh')
        os.system('blockMesh > log.blockMesh')
        print ('Running stitchMesh')
        os.system('chmod 777 stitchMeshSc.sh')
        os.system('bash stitchMeshSc.sh')
        if isOnWallIn0org == 0:
            print ('Making in/onWallCoatFields')
            os.system('python genGeomFPMAll.py')
        print ('Running createPatch')
        os.system('createPatch -overwrite > log.createPatch')
        os.system('paraFoam -touch')
        print ('Creating 0 dir')
        os.system('mkdir 0')
        os.system('cp -rf 0.org/* 0')

    
    # -- just if first run #############################################
    if startTemp == temp:
        # -- set number of iterations for flow in controlDict###########
        print ('Updating controlDict endtime to %d'%(int(nIterFlow)))
        fileHn  = open('system/controlDict','r')
        lineLst = fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'endTime\t%d;\n'%int(nIterFlow)
            if lineLst[ind].find('application') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'application\tporousSimpleFoam;\n'
            if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                spLine = lineLst[ind].split('\t')
                intSpLine = int(spLine[1].replace(';\n',''))
                lineLst[ind] = 'writeInterval \t %d;\n' %int(nIterFlow)
        fileHn.close()
        with open('system/controlDict','w') as fileHn:
            fileHn.writelines(lineLst)
        fileHn.close()
    # --other runs######################################################
    else:
        # --highest results dir#########################################
        dirsI = []
        dirsS = os.listdir('./')
        for Dir in dirsS:
            if isInt(Dir):
                dirsI.append(int(Dir))
        
        # -- set number of iterations for flow in controlDict###########
        print ('Updating controlDict endtime to %d'%(int(nIterFlowThen)+int(max(dirsI))))
        fileHn  = open('system/controlDict','r')
        lineLst = fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'endTime\t%d;\n'%(int(nIterFlowThen)+int(max(dirsI)))
            if lineLst[ind].find('application') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'application\tporousSimpleFoam;\n'
            if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                spLine = lineLst[ind].split('\t')
                intSpLine = int(spLine[1].replace(';\n',''))
                lineLst[ind] = 'writeInterval \t %d;\n' %(int(nIterFlowThen))
        fileHn.close()
        with open('system/controlDict','w') as fileHn:
            fileHn.writelines(lineLst)
        fileHn.close()
    
    # -- read and update inlet temperature in the field ############
    # ~print ('Updating T field to %5.2f'%temp)
    # ~fileHn  = open('%d/T'%max(dirsI),'r')
    # ~lineLst = fileHn.readlines()
    # ~for ind in range(len(lineLst)):
        # ~if lineLst[ind].find('internalField') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~#~ lineLst[ind] =  ' '.join(lineLst[ind].split(' ')[1:len(lineLst[ind].split(' '))]) + ' %5.2f;\n'%temp
            # ~lineLst[ind] = '\tinternalField   uniform %5.3f;\n'%temp
    # ~fileHn.close()
    # ~with open('0.org/T','w') as fileHn:
        # ~fileHn.writelines(lineLst)
    # ~fileHn.close()
    
    # --update and copy p field (p --> kinematic p)#################
    # ~print ('Updating p field in %d'%max(dirsI))
    # ~fileHn  = open('%d/p'%max(dirsI),'r')
    # ~lineLst = fileHn.readlines()
    # ~pNum = False
    # ~for ind in range(len(lineLst)):
        # ~if lineLst[ind].find('dimensions') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = 'dimensions\t[0 2 -2 0 0 0 0];\n'
        # ~if lineLst[ind].find(')\n') >= 0:
            # ~pNum = False
        # ~if pNum:
            # ~lineLst[ind] = '%7.8f\n'%((float(lineLst[ind])-101325.0)/rho)
        # ~if lineLst[ind].find('(\n') >= 0:
            # ~pNum = True
        # ~if lineLst[ind].find('p0') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = '\t\tp0\t\tuniform 0;\n'   
        # ~if lineLst[ind].find('value') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = '\t\tvalue\t\tuniform 0;\n'   
        # ~if lineLst[ind].find('rho') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = '\t\trho\t\tnone;\n'  
    # ~fileHn.close()
    # ~with open('%d/p'%max(dirsI),'w') as fileHn:
        # ~fileHn.writelines(lineLst)
    # ~fileHn.close()

    # --update phi field (divide by rho)############################
    # ~print ('Updating phi field in %d'%max(dirsI))
    # ~fileHn  = open('%d/phi'%max(dirsI),'r')
    # ~lineLst = fileHn.readlines()
    # ~pNum = False
    # ~for ind in range(len(lineLst)):
        # ~if lineLst[ind].find('dimensions') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = 'dimensions\t[0 3 -1 0 0 0 0];\n'
        # ~if lineLst[ind].find(')\n') >= 0:
            # ~pNum = False
        # ~if pNum:
            # ~lineLst[ind] = '%3.10f\n'%float(float(lineLst[ind])/rho)
        # ~if lineLst[ind].find('(\n') >= 0:
            # ~pNum = True
    # ~fileHn.close()
    # ~with open('%d/phi'%max(dirsI),'w') as fileHn:
        # ~fileHn.writelines(lineLst)
    # ~fileHn.close()
    
    # --highest results dir#############################################
    # ~dirsI = []
    # ~dirsS = os.listdir('./')
    # ~for Dir in dirsS:
        # ~if isInt(Dir):
            # ~dirsI.append(int(Dir))
    # ~prevHigh = max(dirsI)
    
    # -- run porousSimpleFoam ##########################################
    # --parallel run not working yet ###################################
    if runParallel:
        os.system('rm -rf log.reconstructPar')
        if not startTemp == temp:
            os.system('cp log.porousSimpleFoam log.porousSimpleFoam%5.2f'%tempLst[ind2-1])
        os.system('rm -rf log.porousSimpleFoam')
        os.system('rm log.renumberMesh')
        os.system('rm log.decomposePar')
        os.system('rm -rf processor*')
        os.system('decomposePar > log.decomposeParDict%5.2f'%tempLst[ind2])
	os.system('./Allrun-parallel')
    else:
        print ('Running porousSimpleFoam with inlet T=%5.2f and nu=%5.8f'%(temp,nu))
        os.system('porousSimpleFoam'' >> log.porousSimpleFoam%5.2f'%temp)
    
    # --highest results dir#############################################
    dirsI = []
    dirsS = os.listdir('./')
    for Dir in dirsS:
        if isInt(Dir):
            dirsI.append(int(Dir))
    
    # -- set number of iterations for kinetics in controlDict###########
    print ('Updating controlDict endtime to %d'%(int(nIterKin)+int(max(dirsI))))
    fileHn  = open('system/controlDict','r')
    lineLst = fileHn.readlines()
    for ind in range(len(lineLst)):
        if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
            lineLst[ind] = 'endTime\t%d;\n'%(int(nIterKin)+int(max(dirsI)))
        if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
            spLine = lineLst[ind].split('\t')
            intSpLine = int(spLine[1].replace(';\n',''))
            lineLst[ind] = 'writeInterval \t %d;\n' %(int(nIterKin-1))
    fileHn.close()
    with open('system/controlDict','w') as fileHn:
        fileHn.writelines(lineLst)
    fileHn.close()
    
    if startTemp == temp:
        # --copy fields from 0.org into highest results directory###########
        print('Copying data from 0.org to %d directory'%max(dirsI))
        os.system('cp 0.org/CO ./%d'%max(dirsI))
        os.system('cp 0.org/O2 ./%d'%max(dirsI))
        os.system('cp 0.org/N2O ./%d'%max(dirsI))
        os.system('cp 0.org/NO ./%d'%max(dirsI))
        # ~os.system('cp 0.org/T ./%d'%max(dirsI))
    else:
        print('Copying data from %d to %d directory'%(prevHigh,max(dirsI)))
        os.system('cp %d/CO ./%d'%(prevHigh,max(dirsI)))
        os.system('cp %d/O2 ./%d'%(prevHigh,max(dirsI)))
        os.system('cp %d/N2O ./%d'%(prevHigh,max(dirsI)))
        os.system('cp %d/NO ./%d'%(prevHigh,max(dirsI)))
        # ~os.system('cp %d/T ./%d'%(prevHigh,max(dirsI)))
    
    os.system('mv -r %d temp%d'%(max(dirsI),temp))
    # --change loadCommonData.H and r_powerLawV3 files that are 
    #   different for each solver
    # ~print ('Updating loadCommonData and powerLaw')
    # ~os.system('cp ./constant/rnTerms/loadCommonDataOM.H ./constant/rnTerms/loadCommonData.H')
    # ~os.system('cp ./constant/rnTerms/r_powerLawV3OM.H ./constant/rnTerms/r_powerLawV3.H')
    
    # --run diffdiffReactConvSysThreeWayMacroSolver#################
    print ('Running diffReactConvSysThreeWayMacro%5.2f'%temp)
    os.system('diffReactConvSysThreeWayMacro''>> log.diffReactConvSysThreeWayMacro%5.2f'%temp)
                
    # --highest results dir#############################################
    prevHigh = max(dirsI)
    dirsI = []
    dirsS = os.listdir('./')
    for Dir in dirsS:
        if isInt(Dir):
            dirsI.append(int(Dir))
            
    # --update and copy p field (kinematic p --> p)#####################
    # ~print ('Updating p field in %d and copy to %d'%(prevHigh,max(dirsI)))
    # ~fileHn  = open('%d/p'%prevHigh,'r')
    # ~lineLst = fileHn.readlines()
    # ~pNum = False
    # ~for ind in range(len(lineLst)):
        # ~if lineLst[ind].find('dimensions') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = 'dimensions\t[1 -1 -2 0 0 0 0];\n'
        # ~if lineLst[ind].find(')\n') >= 0:
            # ~pNum = False
        # ~if pNum:
            # ~lineLst[ind] = '%7.8f\n'%(float(lineLst[ind])*rho+101325)
        # ~if lineLst[ind].find('(\n') >= 0:
            # ~pNum = True
        # ~if lineLst[ind].find('p0') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = '\t\tp0\t\tuniform 101325;\n'   
        # ~if lineLst[ind].find('value') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = '\t\tvalue\t\tuniform 101325;\n'   
        # ~if lineLst[ind].find('rho') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = '\t\trho\t\trho;\n'  
    # ~fileHn.close()
    # ~with open('%d/p'%prevHigh,'w') as fileHn:
        # ~fileHn.writelines(lineLst)
    # ~fileHn.close()
    # ~os.system('cp %d/p %d'%(prevHigh,max(dirsI)))
    
    # --update phi field (multiply by rho)##############################
    # ~print ('Updating phi field in %d'%(max(dirsI)))
    # ~fileHn  = open('%d/phi'%max(dirsI),'r')
    # ~lineLst = fileHn.readlines()
    # ~pNum = False
    # ~for ind in range(len(lineLst)):
        # ~if lineLst[ind].find('dimensions') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = 'dimensions\t[1 0 -1 0 0 0 0];\n'
        # ~if lineLst[ind].find(')\n') >= 0:
            # ~pNum = False
        # ~if pNum:
            # ~lineLst[ind] = '%7.8f\n'%float(float(lineLst[ind])*rho)
        # ~if lineLst[ind].find('(\n') >= 0:
            # ~pNum = True
    # ~fileHn.close()
    # ~with open('%d/phi'%max(dirsI),'w') as fileHn:
        # ~fileHn.writelines(lineLst)
    # ~fileHn.close()
    
    #copy T and other fields ###########################################
    # ~if startTemp == temp:
        # ~print ('Copying T field from 0.org to %d'%max(dirsI))
        # ~os.system('cp 0.org/T ./%d'%max(dirsI))
    os.system('cp %d/p ./%d'%(prevHigh,max(dirsI)))
    # ~else:
        # ~print ('Copying T and p fields from %d to %d'%(prevHigh,max(dirsI)))
        # ~os.system('cp %d/T ./%d'%(prevHigh,max(dirsI)))
        # ~os.system('cp %d/p ./%d'%(prevHigh,max(dirsI)))
            
    # --change loadCommonData.H and r_powerLawV3 files that are 
    #   different for each solver
    # ~print ('Updating loadCommonData and powerLaw')
    # ~os.system('cp ./constant/rnTerms/loadCommonDataNM.H ./constant/rnTerms/loadCommonData.H')
    # ~os.system('cp ./constant/rnTerms/r_powerLawV3NM.H ./constant/rnTerms/r_powerLawV3.H')
    
    # -- set number of iterations for macroMonolith in controlDict######
    # ~print ('Updating controlDict endtime to %d'%(int(nIterAll)+int(max(dirsI))))
    # ~fileHn  = open('system/controlDict','r')
    # ~lineLst = fileHn.readlines()
    # ~for ind in range(len(lineLst)):
        # ~if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
            # ~lineLst[ind] = 'endTime\t%d;\n'%(int(nIterAll)+int(max(dirsI)))
        # ~if lineLst[ind].find('application') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
                # ~lineLst[ind] = 'application\tmacroMonolithFoam;\n'
    # ~fileHn.close()
    # ~with open('system/controlDict','w') as fileHn:
        # ~fileHn.writelines(lineLst)
    # ~fileHn.close()
    
            
    # --run macroMonolithFOAM solver####################################
    # --parallelrun not working yet#####################################
    # ~if runParallel:
        # ~os.system('rm -rf log.reconstructPar')
        # ~os.system('rm log.renumberMesh')
        # ~os.system('rm -rf log.macroMonolithFoam')
        # ~os.system('rm log.decomposePar')
        # ~os.system('rm -rf processor*')
        # ~os.system('./Allrun-parallel')
    # ~else:
        # ~print ('Running macroMonolithFoam%5.2f'%temp)
        # ~os.system('macroMonolithFoam''>> log.macroMonolithFoam%5.2f'%temp)
    
    # --postprocessing #################################################
    os.system("postProcess -func 'patchAverage(name=outlet,CO)' -latestTime > log.patchAverageOutletCO")
    os.system("postProcess -func 'patchAverage(name=outlet,NO)' -latestTime > log.patchAverageOutletNO")
    os.system("postProcess -func 'patchAverage(name=outlet,O2)' -latestTime > log.patchAverageOutletO2")
    os.system("postProcess -func 'patchAverage(name=outlet,N2O)' -latestTime > log.patchAverageOutletN2O")
    os.system("postProcess -func 'patchAverage(name=inlet,p)' -latestTime > log.patchAverageInletP")
    
    # --highest results dir#############################################
    prevHigh = max(dirsI)
    dirsI = []
    dirsS = os.listdir('./')
    for Dir in dirsS:
        if isInt(Dir):
            dirsI.append(int(Dir))
            
    
    # --copy data for next temperature##################################
  #  print ('Copy data into %d => starting directory for next temp'%int(max(dirsI)+1))
  #  os.system('mkdir %d'%int(max(dirsI)+1))
  #  os.system('cp %d/CO ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  os.system('cp %d/O2 ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  os.system('cp %d/N2O ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  os.system('cp %d/NO ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  # ~os.system('cp %d/T ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  os.system('cp %d/U ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  os.system('cp %d/p ./%d'%(int(max(dirsI)),max(dirsI)+1))
  #  os.system('cp %d/phi ./%d'%(int(max(dirsI)),max(dirsI)+1))
    
    # --gather the data into outConcLst including pressure##############
    outConcLst = []
    fileNm = 'log.patchAverageInletP'
    fileHn = open(fileNm,'r')
    lineLst= fileHn.readlines()
    for ind in range(len(lineLst)):
        if lineLst[ind].find('%s'%patchName) >= 0:
            outConcLst.append(float(lineLst[ind].split(' ')[-1]))
    fileHn.close()        
    for specie in speciesLst:
        fileNm = 'log.patchAverageOutlet%s'%(specie)
        fileHn = open(fileNm,'r')
        lineLst= fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('%s'%patchName) >= 0:
                outConcLst.append(1e6*float(lineLst[ind].split(' ')[-1]))
        fileHn.close()
        
    # -- write the results to the final file
    # -- get the size of the output file
    try:
        fileHn = open(outFile,'r')
        startOutFile = False
        fileHn.close()
    except:
        startOutFile = True
    fileHn  = open(outFile,'a')
    if startOutFile:
        fileHn.write('temp\tp\t')
        for specie in speciesLst:
            fileHn.write('%s\t'%specie)
        fileHn.write('\n')

    fileHn.write('%5.2f\t'%temp)
    for outConc in outConcLst:
        fileHn.write('%12.8e\t'%outConc)
    fileHn.write('\n')
    fileHn.close()

endTime = time.time()
print('Total program execution time: %7.2f s\n'%(endTime-startTime))

