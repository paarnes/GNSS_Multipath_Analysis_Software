"""
This module create writes the results to a text file.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import os
import logging
logger = logging.getLogger(__name__)

def writeOutputFile(outputFilename, outputDir, analysisResults, includeResultSummary, includeCompactSummary,\
    includeObservationOverview, includeLLIOverview):

    """
    Function that write out the results of the analysis.

    """

    if outputDir is None or outputDir == "":
        outputDir = 'Outputs_Files'

    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)

    ## - Create full path for storing resultfile
    outputFilename = os.path.join(outputDir, outputFilename)

    ## -- Extracting data
    rinex_obs_filename      = analysisResults['ExtraOutputInfo']['rinex_obs_filename']
    sp3_filename            = analysisResults['ExtraOutputInfo'].get('SP3_filename',None) # added 23.02.2023
    broad_filename          = analysisResults['ExtraOutputInfo'].get('rinex_nav_filename',None) # added 23.02.2023
    markerName              = analysisResults['ExtraOutputInfo']['markerName']
    rinexVersion            = analysisResults['ExtraOutputInfo']['rinexVersion']
    rinexProgr              = analysisResults['ExtraOutputInfo']['rinexProgr']
    recType                 = analysisResults['ExtraOutputInfo']['recType']

    rinex_rec_approx_pos    = analysisResults['ExtraOutputInfo']['Rinex_Receiver_Approx_Pos']
    estimated_approx_pos    = analysisResults['ExtraOutputInfo'].get('Estimated_Receiver_Approx_Pos', None)
    esti_approx_pos_stats   = analysisResults['ExtraOutputInfo'].get('Estimated_Receiver_Approx_Pos_stats', None)

    tFirstObs               = analysisResults['ExtraOutputInfo']['tFirstObs']
    tLastObs                = analysisResults['ExtraOutputInfo']['tLastObs']
    tInterval               = analysisResults['ExtraOutputInfo']['tInterval']
    GLO_Slot2ChannelMap     = analysisResults['ExtraOutputInfo']['GLO_Slot2ChannelMap']
    nClockJumps             = analysisResults['ExtraOutputInfo']['nClockJumps']
    stdClockJumpInterval    = analysisResults['ExtraOutputInfo']['stdClockJumpInterval']
    meanClockJumpInterval   = analysisResults['ExtraOutputInfo']['meanClockJumpInterval']
    ionLimit                = analysisResults['ExtraOutputInfo']['ionLimit']
    phaseCodeLimit          = analysisResults['ExtraOutputInfo']['phaseCodeLimit']
    elevation_cutoff        = analysisResults['ExtraOutputInfo']['elevation_cutoff'] # added 23.02.2023



    ## Extract systems in current analysis
    GNSSsystems = analysisResults['GNSSsystems']
    GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C']))
    GNSS_Name2Code = {key:val for key, val in GNSS_Name2Code.items() if val in GNSSsystems} # only the systems for current analysis
    ## Replace letter with whole system name
    if 'G' in GNSSsystems:
        GNSSsystems[GNSSsystems.index('G')] = 'GPS'
    if 'R' in GNSSsystems:
        GNSSsystems[GNSSsystems.index('R')] = 'GLONASS'
    if 'E' in GNSSsystems:
        GNSSsystems[GNSSsystems.index('E')] = 'Galileo'
    if 'C' in GNSSsystems:
        GNSSsystems[GNSSsystems.index('C')] = 'BeiDou'


    nGNSSsystems = len(GNSSsystems)
    YesNoMap = {1 : 'Yes',
                0 : 'No'}


    GPSBandNameMap      = dict(zip([1, 2, 5], ['L1', 'L2', 'L5']))
    GLONASSBandNameMap  = dict(zip([1, 2, 3, 4, 6], ['G1', 'G2', 'G3', 'G1a', 'G2a']))
    GalileoBandNameMap  = dict(zip([1, 5, 6, 7, 8], ['E1', 'E5a', 'E6', 'E5b', 'E5(a+b)']))
    BeiDouBandNameMap   = dict(zip([1, 2, 5, 6, 7, 8], ['B1', 'B1-2', 'B2a', 'B3', 'B2b', 'B2(a+b)']))
    GNSSBandNameMap     = dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], [GPSBandNameMap, GLONASSBandNameMap, GalileoBandNameMap, BeiDouBandNameMap]))

    GPSBandFreqMap      = dict(zip([1, 2, 5], ['1575.42', '1227.60', '1176.45']))
    GLONASSBandFreqMap  = dict(zip([1, 2, 3, 4, 6], ['1602 + k*9/16', '1246 + k*7/16', '1202.025', '1600.995', '1248.06']))
    GalileoBandFreqMap  = dict(zip([1, 5, 6, 7, 8], ['1575.42', '1176.45', '1278.75', '1207.140', '1191.795']))
    BeiDouBandFreqMap   = dict(zip([1, 2, 5, 6, 7, 8], ['1575.42', '1561.098', '1176.45', '1268.52', '1207.140', '1191.795']))
    GNSSBandFreqMap     = dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], [GPSBandFreqMap, GLONASSBandFreqMap, GalileoBandFreqMap, BeiDouBandFreqMap]))

    ## -- Check if any LLI indicators at all
    LLI_Active = 0
    for i in range(0,nGNSSsystems):
        current_sys_struct = analysisResults[list(GNSS_Name2Code.keys())[i]]
        nBands = current_sys_struct['nBands']
        for j in range(0,nBands):
            current_band_struct = current_sys_struct[current_sys_struct['Bands'][j]]
            nCodes = current_band_struct['nCodes']
            for k in range(0,nCodes):
                try:
                    current_code_struct = current_band_struct[current_band_struct['Codes'][k]]
                except:
                    continue # If noe code available

                if current_code_struct['LLI_slip_distribution']['n_slips_Tot'] > 0:
                    LLI_Active = 1
                    continue
            if LLI_Active:
                continue
        if LLI_Active:
            continue
    if not LLI_Active:
        includeLLIOverview = 0

    ## -- HEADER
    fid = open(outputFilename, 'w+')

    fid.write('GNSS_MultipathAnalysis\n')
    fid.write('Software version: 1.5.0\n')
    fid.write('Last software version release: 08/01/2025\n\n')
    fid.write('Software developed by Per Helge Aarnes (per.helge.aarnes@gmail.com) \n\n')
    fid.write('RINEX observation filename:\t\t %s\n' % (rinex_obs_filename))
    if sp3_filename is not None:
        fid.write('SP3 filename:\t\t\t\t\t %s\n' % (','.join(sp3_filename)))
    else:
        fid.write('Broadcast navigation filename:\t %s\n' % (','.join(broad_filename)))
    fid.write('RINEX version:\t\t\t\t\t %s\n' % (rinexVersion.strip()))
    fid.write('RINEX converting program:\t\t %s\n' % (rinexProgr))
    fid.write('Rec. approx position (RINEX):    %s\n' % (rinex_rec_approx_pos))

    if estimated_approx_pos:
        std_est_pos = (
            esti_approx_pos_stats["Standard Deviations"]["Sx"],
            esti_approx_pos_stats["Standard Deviations"]["Sy"],
            esti_approx_pos_stats["Standard Deviations"]["Sz"],
        )
        fid.write('Est. approx position: \t\t\t %s\n' % str(estimated_approx_pos))
        fid.write('St.dev of the est. position: \t %s\n' % ', '.join(map(str, std_est_pos)))

    # Convert to 1D and extract scalars
    tFirstObs_flat = tFirstObs.flatten()
    tLastObs_flat = tLastObs.flatten()

    fid.write('Marker name:\t\t\t\t\t %s\n' % (markerName))
    fid.write('Receiver type:\t\t\t\t\t %s\n' % (recType))
    fid.write('Date of observation start:\t\t %4d/%d/%d %d:%d:%.2f \n' % (tFirstObs_flat[0],tFirstObs_flat[1],tFirstObs_flat[2],tFirstObs_flat[3],tFirstObs_flat[4],tFirstObs_flat[5]))
    fid.write('Date of observation end:\t\t %4d/%d/%d %d:%d:%.2f \n'   % (tLastObs_flat[0],tLastObs_flat[1],tLastObs_flat[2],tLastObs_flat[3],tLastObs_flat[4],tLastObs_flat[5]))
    fid.write('Observation interval [seconds]:\t %d\n' % (tInterval))
    fid.write('Elevation angle cutoff [degree]: %d\n' % (elevation_cutoff))
    fid.write('Number of receiver clock jumps:\t %d\n' % (nClockJumps))
    fid.write('Average clock jumps interval:\t %s (std: %.2f seconds)\n\n' % (str(meanClockJumpInterval), stdClockJumpInterval))

    fid.write('Critical cycle slip limits [m/s]:\n')
    fid.write('- Ionospheric delay:\t\t\t%6.3f\n'% (ionLimit))
    fid.write('- Phase-code combination:\t\t%6.3f\n\n' % (phaseCodeLimit))
    fid.write('GNSS systems presents in RINEX observation file:\n')


    for i in range(0,nGNSSsystems):
        fid.write('- %s\n' % (analysisResults['GNSSsystems'][i]))


    if not LLI_Active:
        fid.write('\n\nNOTE: As there were no "Loss-of-Lock" indicators in RINEX observation file,\n. No information concerining "Loss-of-Lock" indicators is included in output file')


    fid.write('\n\nUser-specified contend included in output file\n')
    fid.write('- Include overview of observations for each satellite:\t\t\t%s\n' % (YesNoMap[includeObservationOverview]))
    fid.write('- Include compact summary of analysis estimates:\t\t\t\t%s\n' % (YesNoMap[includeCompactSummary]))
    fid.write('- Include detailed summary of analysis\n   estimates, including for each individual satellite:\t\t\t%s\n' % (YesNoMap[includeResultSummary]))
    fid.write('- Include information about "Loss-of-Lock"\n   indicators in detailed summary:\t\t\t\t\t\t\t\t%s\n' % (YesNoMap[includeLLIOverview]))

    fid.write('\n\n======================================================================================================================================================================================================================================================================================================================================================\n')
    fid.write('======================================================================================================================================================================================================================================================================================================================================================\n\n')
    fid.write('END OF HEADER\n\n\n')


    ## -- COMPLETENESS OVERVIEW
    if includeObservationOverview:
        fid.write( '\n\n\n\n======================================================================================================================================================================================================================================================================================================================================================')
        fid.write( '\n======================================================================================================================================================================================================================================================================================================================================================\n\n')
        fid.write( 'OBSERVATION COMPLETENESS OVERVIEW\n\n\n')
        for i in range(0,nGNSSsystems):
            if GNSSsystems[i] == 'GPS':
                fid.write( 'GPS Observation overview\n')
                nSat = len(analysisResults['GPS']['observationOverview'])
                fid.write( ' ___________________________________________________________________________________________________\n')
                fid.write( '|  PRN   |        L1 Observations          |             L2 Observations          | L5 Observations |\n')
                for PRN in range(0,nSat):
                    PRN = PRN +1

                    if not all([analysisResults['GPS']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                              analysisResults['GPS']['observationOverview']['Sat_' + str(PRN)]['Band_2'],\
                              analysisResults['GPS']['observationOverview']['Sat_' + str(PRN)]['Band_5']]) == "":

                        fid.write(  '|________|_________________________________|______________________________________|_________________|\n')
                        fid.write( '|%8s|%33s|%38s|%17s|\n' % ( \
                            GNSS_Name2Code[analysisResults['GNSSsystems'][i]] + str(PRN), \
                            analysisResults['GPS']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                            analysisResults['GPS']['observationOverview']['Sat_' + str(PRN)]['Band_2'],\
                            analysisResults['GPS']['observationOverview']['Sat_' + str(PRN)]['Band_5']))

                fid.write( '|________|_________________________________|______________________________________|_________________|\n\n\n')

            elif GNSSsystems[i] == 'GLONASS':

                fid.write(  'GLONASS Observation overview\n')
                nSat =  len(analysisResults['GLONASS']['observationOverview'])
                fid.write(  ' ________________________________________________________________________________________________________________________\n')
                fid.write(  '| Sat ID | Frequency Channel | G1 Observations | G2 Observations | G3 Observations | G1a Observations | G2a Observations |\n')
                for PRN in list(GLO_Slot2ChannelMap.keys()):
                    if not all([analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                           analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_2'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_3'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_4'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_6']]) == "":

                        fid.write(  '|________|___________________|_________________|_________________|_________________|__________________|__________________|\n')
                        fid.write(  '|%8s|%19d|%17s|%17s|%17s|%18s|%18s|\n' % (\
                            GNSS_Name2Code[analysisResults['GNSSsystems'][i]]+ str(PRN), \
                            GLO_Slot2ChannelMap[PRN],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_2'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_3'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_4'],\
                            analysisResults['GLONASS']['observationOverview']['Sat_' + str(PRN)]['Band_1']))

                fid.write(  '|________|___________________|_________________|_________________|_________________|__________________|__________________|\n\n\n')

            elif GNSSsystems[i] == 'Galileo':

                fid.write(  'Galileo Observation overview\n')
                nSat = len(analysisResults['Galileo']['observationOverview'].keys())
                fid.write(  ' _________________________________________________________________________________________________________\n')
                fid.write(  '|  PRN   | E1 Observations | E5a Observations | E6 Observations | E5b Observations | G5(a+b) Observations |\n')
                for PRN in range(0,nSat):
                    PRN = PRN + 1
                    if not all([analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                            analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_5'],\
                            analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_6'],\
                            analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_7'],\
                            analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_8']]) == "":

                        fid.write(   '|________|_________________|__________________|_________________|__________________|______________________|\n')
                        fid.write(   '|%8s|%17s|%18s|%17s|%18s|%22s|\n' % (\
                            GNSS_Name2Code[analysisResults['GNSSsystems'][i]]+ str(PRN), \
                            analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                             analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_5'],\
                             analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_6'],\
                             analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_7'],\
                             analysisResults['Galileo']['observationOverview']['Sat_' + str(PRN)]['Band_8']))

                fid.write( '|________|_________________|__________________|_________________|__________________|______________________|\n\n\n')

            elif GNSSsystems[i] =='BeiDou':

                fid.write(  'BeiDou Observation overview\n')
                nSat = len(analysisResults['BeiDou']['observationOverview'].keys())
                fid.write(  ' ______________________________________________________________________________________________________________________________\n')
                fid.write(  '|  PRN   | B1 Observations | E1-2 Observations | B2a Observations | B3 Observations  | B2b Observations | B2(a+b) Observations |\n')
                for PRN in range(0,nSat):
                    PRN = PRN + 1
                    if not all([analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_2'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_5'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_6'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_7'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_8']]) == "":

                        fid.write(  '|________|_________________|___________________|__________________|__________________|__________________|______________________|\n')
                        fid.write( '|%8s|%17s|%19s|%18s|%18s|%18s|%22s|\n' % (\
                            GNSS_Name2Code[analysisResults['GNSSsystems'][i]]+ str(PRN), \
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_1'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_2'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_5'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_6'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_7'],\
                            analysisResults['BeiDou']['observationOverview']['Sat_' + str(PRN)]['Band_8']))

                fid.write( '|________|_________________|___________________|__________________|__________________|__________________|______________________|\n\n\n')

        fid.write(  '======================================================================================================================================================================================================================================================================================================================================================\n')
        fid.write( '======================================================================================================================================================================================================================================================================================================================================================\n')
        fid.write( 'END OF OBSERVATION COMPLETENESS OVERVIEW\n\n\n\n\n')

    ## -- Compact Code analysis summary

    if includeCompactSummary:
        fid.write( '\n\n\n\n======================================================================================================================================================================================================================================================================================================================================================')
        fid.write(  '\n======================================================================================================================================================================================================================================================================================================================================================\n\n')
        fid.write(  'ANALYSIS RESULTS SUMMARY (COMPACT)\n\n\n')


        for i in range(0,nGNSSsystems):
            curr_sys = GNSSsystems[i]
            current_sys_struct = analysisResults[analysisResults['GNSSsystems'][i]]
            nBands_current_sys = current_sys_struct['nBands']
            current_BandFreqMap = GNSSBandFreqMap[GNSSsystems[i]]
            current_BandNameMap = GNSSBandNameMap[GNSSsystems[i]]

            headermsg                       = '|                                             |'
            rmsMultiPathmsg                 = '|RMS multipath[meters]                        |'
            rmsMultiPathmsg_weighted        = '|Weighted RMS multipath[meters]               |'
            nSlipsmsg                       = '|N ambiguity slips periods                    |'
            slipRatiomsg                    = '|Ratio of N slip periods/N obs epochs [%]     |'
            nSlipsOver10msg                 = '|N slip periods, elevation angle > 10 degrees |'
            nSlipsUnder10msg                = '|N slip periods, elevation angle < 10 degrees |'
            nSlipsNaNmsg                    = '|N slip periods, elevation angle not computed |'
            topline                         = ' _____________________________________________'
            bottomline                      = '|_____________________________________________|'

            fid.write(  '\n\n\n\n')
            fid.write(  '%s ANALYSIS SUMMARY\n\n' % (GNSSsystems[i]))
            for j in range(0,nBands_current_sys):
                bandName = current_sys_struct['Bands'][j]
                current_band_struct = current_sys_struct[bandName]
                nCodes_current_band = current_band_struct['nCodes']
                for k in range(0,nCodes_current_band):
                    codeName = current_band_struct['Codes'][k]
                    try:
                        current_code_struct = current_band_struct[codeName]
                    except:
                        logger.warning(f"INFO(GNSS_MultipathAnalysis): No estimates to put in report for {codeName} for {curr_sys}")
                        continue

                    topline                     = topline + '_________'
                    bottomline                  = bottomline + '________|'
                    headermsg                   = headermsg + '%8s|' % (codeName)
                    rmsMultiPathmsg             = rmsMultiPathmsg + '%8.3f|' % (current_code_struct['rms_multipath_range1_averaged'])
                    rmsMultiPathmsg_weighted    = rmsMultiPathmsg_weighted + '%8.3f|' % (current_code_struct['elevation_weighted_average_rms_multipath_range1'])
                    slipRatiomsg                = slipRatiomsg +  '%8.3f|' % (100*current_code_struct['range1_slip_distribution']['n_slips_Tot']/current_code_struct['nRange1Obs'])
                    nSlipsmsg                   = nSlipsmsg + '%8d|' % (current_code_struct['range1_slip_distribution']['n_slips_Tot'])
                    nSlipsOver10msg             = nSlipsOver10msg + '%8d|' % ( \
                        sum([current_code_struct['range1_slip_distribution']['n_slips_10_20'], current_code_struct['range1_slip_distribution']['n_slips_20_30'], \
                             current_code_struct['range1_slip_distribution']['n_slips_30_40'], current_code_struct['range1_slip_distribution']['n_slips_40_50'], \
                             current_code_struct['range1_slip_distribution']['n_slips_over50']]))

                    nSlipsUnder10msg            = nSlipsUnder10msg + '%8d|' % (current_code_struct['range1_slip_distribution']['n_slips_0_10'])
                    nSlipsNaNmsg                = nSlipsNaNmsg + '%8d|' % (current_code_struct['range1_slip_distribution']['n_slips_NaN'])

            fid.write(topline + '\n')
            fid.write(headermsg + '\n')
            fid.write(bottomline +'\n')
            fid.write(rmsMultiPathmsg + '\n')
            fid.write(bottomline + '\n')
            fid.write(rmsMultiPathmsg_weighted + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlipsmsg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlipsOver10msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlipsUnder10msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlipsNaNmsg + '\n')
            fid.write(bottomline + '\n')
            fid.write(slipRatiomsg + '\n')
            fid.write(bottomline + '\n')



            ## -- Summary for cycle slip detected with both ionospheric residuals and code-phase difference
            fid.write(  '\n\n') # make some space
            headermsg                       = '|                                             |'
            nSlipsmsg                       = '|N detected cycle slips                       |'
            slipRatiomsg                    = '|Ratio of N cycle slips/N obs epochs [%]      |'
            nSlipsUnder10msg                = '|N cycle slip, elevation angle < 10 degrees   |'
            nSlips10_20msg                  = '|N cycle slip, elevation angle 10-20 degrees  |'
            nSlips20_30msg                  = '|N cycle slip, elevation angle 20-30 degrees  |'
            nSlips30_40msg                  = '|N cycle slip, elevation angle 30-40 degrees  |'
            nSlips40_50msg                  = '|N cycle slip, elevation angle 40-50 degrees  |'
            nSlipsOver50msg                 = '|N cycle slip, elevation angle > 50 degrees   |'
            nSlipsNaNmsg                    = '|N cycle slip, elevation angle not computed   |'
            topline                         = ' _____________________________________________'
            bottomline                      = '|_____________________________________________|'

            fid.write('\n')
            fid.write(  '%s: DETECTED CYCLE SLIPS IN TOTAL FOR THE SIGNAL COMBINATION (IONOSPHERIC RESIDUALS & CODE-PHASE COMBINATION)\n' % (GNSSsystems[i]))
            for j in range(0,nBands_current_sys):
                bandName = current_sys_struct['Bands'][j]
                current_band_struct = current_sys_struct[bandName]
                nCodes_current_band = current_band_struct['nCodes']
                for k in range(0,nCodes_current_band):
                    codeName = current_band_struct['Codes'][k]
                    try:
                        current_code_struct = current_band_struct[codeName]
                    except:
                        continue

                    topline                     = topline + '_________'
                    bottomline                  = bottomline + '________|'
                    headermsg                   = headermsg + '%8s|' % (codeName)
                    slipRatiomsg                = slipRatiomsg +  '%8.3f|' % (100*current_code_struct['cycle_slip_distribution']['n_slips_Tot']/current_code_struct['nRange1Obs'])
                    nSlipsmsg                   = nSlipsmsg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_Tot'])
                    nSlipsUnder10msg            = nSlipsUnder10msg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_0_10'])
                    nSlips10_20msg              = nSlips10_20msg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_10_20'])
                    nSlips20_30msg              = nSlips20_30msg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_20_30'])
                    nSlips30_40msg              = nSlips30_40msg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_30_40'])
                    nSlips40_50msg              = nSlips40_50msg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_40_50'])
                    nSlipsOver50msg             = nSlipsOver50msg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_over50'])
                    # nSlipsNaNmsg                = nSlipsNaNmsg + '%8d|' % (current_code_struct['cycle_slip_distribution']['n_slips_NaN'])

            fid.write(topline + '\n')
            fid.write(headermsg + '\n')
            fid.write(bottomline +'\n')
            fid.write(nSlipsmsg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlipsUnder10msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlips10_20msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlips20_30msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlips30_40msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlips40_50msg + '\n')
            fid.write(bottomline + '\n')
            fid.write(nSlipsOver50msg + '\n')
            fid.write(bottomline + '\n')
            # fid.write(nSlipsNaNmsg + '\n')
            # fid.write(bottomline + '\n')
            fid.write(slipRatiomsg + '\n')
            fid.write(bottomline + '\n')


        fid.write(   '\n======================================================================================================================================================================================================================================================================================================================================================\n')
        fid.write(   'END OF ANALYSIS RESULTS SUMMARY (COMPACT)\n\n\n\n\n')



    ## -- Code analysis
    if includeResultSummary:
        for i in range(0,nGNSSsystems):
            current_sys = GNSSsystems[i]
            current_sys_struct = analysisResults[GNSSsystems[i]]
            nBands_current_sys = current_sys_struct['nBands']
            current_BandFreqMap = GNSSBandFreqMap[GNSSsystems[i]]
            current_BandNameMap = GNSSBandNameMap[GNSSsystems[i]]

            fid.write('\n\n\n\n======================================================================================================================================================================================================================================================================================================================================================')
            fid.write('\n======================================================================================================================================================================================================================================================================================================================================================\n\n')
            fid.write('BEGINNING OF %s ANALYSIS\n\n' % (analysisResults['GNSSsystems'][i]))
            fid.write('Amount of carrier bands analysed: %d \n' % nBands_current_sys)
            for j in range(0,nBands_current_sys):
                bandName = current_sys_struct['Bands'][j]
                current_band_struct = current_sys_struct[bandName]
                nCodes_current_band = current_band_struct['nCodes']
                fid.write( '\n\n======================================================================================================================================================================================================================================================================================================================================================\n\n')
                if int(bandName[-1]) in current_BandNameMap.keys():
                    fid.write( '%s (%s)\n\n' % (analysisResults[GNSSsystems[i]]['Bands'][j], current_BandNameMap[int(bandName[-1])]))
                else:
                    continue
                fid.write( 'Frequency of carrier band [MHz]:\t\t\t\t\t %s\n' % (current_BandFreqMap[int(bandName[-1])]))
                fid.write( 'Amount of code signals analysed in current band:\t %d \n' % (nCodes_current_band))
                for k in range(0,nCodes_current_band):
                    try:
                        current_code_struct = current_band_struct[current_band_struct['Codes'][k]]
                    except:
                        break
                    fid.write( '\n------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n')
                    fid.write(  'Code signal:\t\t\t\t\t\t\t\t\t %s\n\n' % current_code_struct['range1_Code'])
                    fid.write(  'Second code signal\n(Utilized for linear combinations):\t\t\t\t %s\n' % (current_code_struct['range2_Code']))
                    fid.write( 'RMS multipath (All SVs) [meters]:\t\t\t\t%6.3f \n' % (current_code_struct['rms_multipath_range1_averaged']))
                    fid.write(  'Weighted RMS multipath (All SVs) [meters]:\t\t%6.3f \n' %  (current_code_struct['elevation_weighted_average_rms_multipath_range1']))
                    fid.write( 'Number of %s observation epochs:\t\t\t\t %d \n' % (current_code_struct['range1_Code'], current_code_struct['nRange1Obs']))
                    fid.write(  'N epochs with multipath estimates:\t\t\t\t %d \n' % (current_code_struct['nEstimates']))
                    fid.write(  'N ambiguity slips on %s signal:\t\t\t\t %d \n'  % (\
                        current_code_struct['range1_Code'], current_code_struct['range1_slip_distribution']['n_slips_Tot']))
                    fid.write(  'Ratio of N slip periods/N %s obs epochs [%%]:\t %.3f\n' % (\
                        current_code_struct['range1_Code'], 100*current_code_struct['range1_slip_distribution']['n_slips_Tot']/current_code_struct['nRange1Obs']))


                    nSat = len(current_code_struct['range1_slip_distribution_per_sat'])

                    if includeLLIOverview:

                        if not current_sys == 'GLONASS':
                            fid.write(  '\nSatellite Overview\n')
                            fid.write( ' _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________\n')
                            fid.write( '|   |    n %s   | n Epochs with |   RMS   | Weighted RMS |  Average Sat. |                           |     Slip Periods/Obs      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |\n' % current_code_struct['range1_Code'])
                            fid.write( '|PRN|Observations|   Multipath   |Multipath|  Multipath   |Elevation Angle|       n Slip Periods      |         Ratio             |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |\n')
                            fid.write( '|   |            |   Estimates   |[meters] |   [meters]   |   [degrees]   |                           |          [%]              |       0-10 degrees        |        10-20 degrees      |        20-30 degrees      |        30-40 degrees      |        40-50 degrees      |        >50 degrees        |        NaN degrees        |\n')
                            fid.write( '|   |            |               |         |              |               |___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|\n')
                            fid.write( '|   |            |               |         |              |               | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  |\n')

                        for PRN in range(0, nSat):
                            if current_code_struct['nEstimates_per_sat'][PRN] > 0:
                                fid.write('|___|____________|_______________|_________|______________|_______________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|\n')
                                fid.write('|%3s|%12d|%15d|%9.3f|%14.3f|%15.3f|%10d|%7d|%8d|%10.3f|%7.3f|%8.3f|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|\n' % (
                                    GNSS_Name2Code[GNSSsystems[i]] + str(PRN),
                                    int(current_code_struct['n_range1_obs_per_sat'][:, PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:, PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:, PRN]),
                                    int(current_code_struct['nEstimates_per_sat'][PRN].item()) if hasattr(current_code_struct['nEstimates_per_sat'][PRN], "item") else int(current_code_struct['nEstimates_per_sat'][PRN]),
                                    float(current_code_struct['rms_multipath_range1_satellitewise'][PRN]),
                                    float(current_code_struct['elevation_weighted_rms_multipath_range1_satellitewise'][PRN]),
                                    float(current_code_struct['mean_sat_elevation_angles'][PRN]),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_Tot']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_Tot']),
                                    100 * int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot']) / (int(current_code_struct['n_range1_obs_per_sat'][:, PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:, PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:, PRN])),
                                    100 * int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_Tot']) / (int(current_code_struct['n_range1_obs_per_sat'][:, PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:, PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:, PRN])),
                                    100 * int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_Tot']) / (int(current_code_struct['n_range1_obs_per_sat'][:, PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:, PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:, PRN])),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_0_10']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_0_10']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_0_10']),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_10_20']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_10_20']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_10_20']),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_20_30']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_20_30']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_20_30']),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_30_40']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_30_40']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_30_40']),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_40_50']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_40_50']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_40_50']),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_over50']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_over50']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_over50']),
                                    int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_NaN']),
                                    int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_NaN']),
                                    int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_NaN'])
                                ))

                            fid.write(  '|___|____________|_______________|_________|______________|_______________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|\n')
                        else:
                            fid.write(  '\nSatellite Overview\n')
                            fid.write(  ' ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________\n')
                            fid.write(  '|      | Frequency |    n %s   | n Epochs with |   RMS   | Weighted RMS |  Average Sat. |                           |        Slip/Obs           |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |       n Slip Periods      |\n' % current_code_struct['range1_Code'])
                            fid.write(  '|Sat ID|  Channel  |Observations|   Multipath   |Multipath|  Multipath   |Elevation Angle|       n Slip Periods      |         Ratio             |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |      Elevation Angle      |\n')
                            fid.write(  '|      |           |            |   Estimates   |[meters] |   [meters]   |   [degrees]   |                           |          [%]              |       0-10 degrees        |        10-20 degrees      |        20-30 degrees      |        30-40 degrees      |        40-50 degrees      |        >50 degrees        |        NaN degrees        |\n')
                            fid.write(  '|      |           |            |               |         |              |               |___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|___________________________|\n')
                            fid.write(  '|      |           |            |               |         |              |               | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  | Analysed |  LLI  |  Both  |\n')
                            # for PRN in range(0,nSat):
                            for PRN in list(GLO_Slot2ChannelMap.keys()):
                                if current_code_struct['nEstimates_per_sat'][PRN] > 0:
                                    fid.write('|______|___________|____________|_______________|_________|______________|_______________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|\n')
                                    fid.write('|%6s|%11d|%12d|%15d|%9.3f|%14.3f|%15.3f|%10d|%7d|%8d|%10.3f|%7.3f|%8.3f|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|%10d|%7d|%8d|\n' % (
                                        GNSS_Name2Code[GNSSsystems[i]] + str(PRN),
                                        int(GLO_Slot2ChannelMap[PRN].item()) if hasattr(GLO_Slot2ChannelMap[PRN], "item") else int(GLO_Slot2ChannelMap[PRN]),
                                        int(current_code_struct['n_range1_obs_per_sat'][:,PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:,PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:,PRN]),
                                        int(current_code_struct['nEstimates_per_sat'][PRN].item()) if hasattr(current_code_struct['nEstimates_per_sat'][PRN], "item") else int(current_code_struct['nEstimates_per_sat'][PRN]),
                                        float(current_code_struct['rms_multipath_range1_satellitewise'][PRN]),
                                        float(current_code_struct['elevation_weighted_rms_multipath_range1_satellitewise'][PRN]),
                                        float(current_code_struct['mean_sat_elevation_angles'][PRN]),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_Tot']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_Tot']),
                                        100 * int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot']) / (int(current_code_struct['n_range1_obs_per_sat'][:,PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:,PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:,PRN])),
                                        100 * int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_Tot']) / (int(current_code_struct['n_range1_obs_per_sat'][:,PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:,PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:,PRN])),
                                        100 * int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_Tot']) / (int(current_code_struct['n_range1_obs_per_sat'][:,PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:,PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:,PRN])),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_0_10']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_0_10']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_0_10']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_10_20']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_10_20']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_10_20']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_20_30']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_20_30']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_20_30']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_30_40']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_30_40']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_30_40']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_40_50']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_40_50']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_40_50']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_over50']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_over50']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_over50']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_NaN']),
                                        int(current_code_struct['LLI_slip_distribution_per_sat'][PRN]['n_slips_NaN']),
                                        int(current_code_struct['slip_distribution_per_sat_LLI_fusion'][PRN]['n_slips_NaN'])
                                    ))

                            fid.write(  '|______|___________|____________|_______________|_________|______________|_______________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|__________|_______|________|\n')

                    else:
                        if not current_sys == 'GLONASS':
                            fid.write(   '\nSatellite Overview\n')
                            fid.write(   ' __________________________________________________________________________________________________________________________________________________________________________________________________________________________________ \n')
                            fid.write(   '|   |    n %s   | n Epochs with |   RMS   | Weighted RMS |  Average Sat. |               | Slip/Obs | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  |\n' % (current_code_struct['range1_Code']))
                            fid.write(   '|PRN|Observations|   Multipath   |Multipath|  Multipath   |Elevation Angle|    n Slip     |  Ratio   | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle |\n')
                            fid.write(   '|   |            |   Estimates   |[meters] |   [meters]   |   [degrees]   |    Periods    |   [%]    |  0-10 degrees   |  10-20 degrees  |  20-30 degrees  |  30-40 degrees  |  40-50 degrees  |   >50 degrees   |   NaN degrees   |\n')

                            for PRN in range(0, nSat):
                                if current_code_struct['nEstimates_per_sat'][PRN] > 0:  # added 21.01.2023 to prevent sat with only nan in resultfile
                                    fid.write('|___|____________|_______________|_________|______________|_______________|_______________|__________|_________________|_________________|_________________|_________________|_________________|_________________|_________________|\n')
                                    n_obs = int(current_code_struct['n_range1_obs_per_sat'][:, PRN].item()) if hasattr(current_code_struct['n_range1_obs_per_sat'][:, PRN], "item") else int(current_code_struct['n_range1_obs_per_sat'][:, PRN])
                                    n_est = int(current_code_struct['nEstimates_per_sat'][PRN].item()) if hasattr(current_code_struct['nEstimates_per_sat'][PRN], "item") else int(current_code_struct['nEstimates_per_sat'][PRN])
                                    slip_tot = int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot'])
                                    slip_ratio = 100 * slip_tot / n_obs if n_obs != 0 else float('nan')
                                    fid.write('|%3s|%12d|%15d|%9.3f|%14.3f|%15.3f|%15d|%10.3f|%17d|%17d|%17d|%17d|%17d|%17d|%17d|\n' % (
                                        GNSS_Name2Code[GNSSsystems[i]] + str(PRN),
                                        n_obs,
                                        n_est,
                                        float(current_code_struct['rms_multipath_range1_satellitewise'][PRN]),
                                        float(current_code_struct['elevation_weighted_rms_multipath_range1_satellitewise'][PRN]),
                                        float(current_code_struct['mean_sat_elevation_angles'][PRN]),
                                        slip_tot,
                                        slip_ratio,
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_0_10']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_10_20']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_20_30']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_30_40']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_40_50']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_over50']),
                                        int(current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_NaN'])
                                    ))

                            fid.write(  '|___|____________|_______________|_________|______________|_______________|_______________|__________|_________________|_________________|_________________|_________________|_________________|_________________|_________________|\n')
                        else:
                            fid.write(  '\nSatellite Overview\n')
                            fid.write(  ' _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________ \n')
                            fid.write(  '|      | Frequency |    n %s   | n Epochs with |   RMS   | Weighted RMS |  Average Sat. |               | Slip/Obs | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  | n Slip Periods  |\n' % current_code_struct['range1_Code'])
                            fid.write(  '|Sat ID|  Channel  |Observations|   Multipath   |Multipath|  Multipath   |Elevation Angle|    n Slip     |  Ratio   | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle | Elevation Angle |\n')
                            fid.write(  '|      |           |            |   Estimates   |[meters] |   [meters]   |   [degrees]   |    Periods    |   [%]    |  0-10 degrees   |  10-20 degrees  |  20-30 degrees  |  30-40 degrees  |  40-50 degrees  |   >50 degrees   |   NaN degrees   |\n')

                            for PRN in list(GLO_Slot2ChannelMap.keys()):
                                if current_code_struct['nEstimates_per_sat'][PRN] > 0: ##added 21.01.2023 to prevent sat with only nan in resultfile
                                    fid.write(   '|______|___________|____________|_______________|_________|______________|_______________|_______________|__________|_________________|_________________|_________________|_________________|_________________|_________________|_________________|\n')
                                    fid.write(   '|%6s|%11d|%12d|%15d|%9.3f|%14.3f|%15.3f|%15d|%10.3f|%17d|%17d|%17d|%17d|%17d|%17d|%17d|\n' % (\
                                        # GNSS_Name2Code[analysisResults[GNSSsystems[i]]] + str(PRN),\
                                        GNSS_Name2Code[GNSSsystems[i]] + str(PRN),\
                                        GLO_Slot2ChannelMap[PRN],\
                                        current_code_struct['n_range1_obs_per_sat'][:,PRN],\
                                        current_code_struct['nEstimates_per_sat'][PRN],\
                                        current_code_struct['rms_multipath_range1_satellitewise'][PRN],\
                                        current_code_struct['elevation_weighted_rms_multipath_range1_satellitewise'][PRN],\
                                        current_code_struct['mean_sat_elevation_angles'][PRN],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot'],\
                                        100*current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_Tot']/current_code_struct['n_range1_obs_per_sat'][:,PRN],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_0_10'],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_10_20'],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_20_30'],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_30_40'],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_40_50'],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_over50'],\
                                        current_code_struct['range1_slip_distribution_per_sat'][PRN]['n_slips_NaN']))

                            fid.write(   '|______|___________|____________|_______________|_________|______________|_______________|_______________|__________|_________________|_________________|_________________|_________________|_________________|_________________|_________________|\n')

            fid.write(   '\n======================================================================================================================================================================================================================================================================================================================================================\n')
            fid.write(   '======================================================================================================================================================================================================================================================================================================================================================\n')
            fid.write(   'END OF %s ANALYSIS\n\n\n\n' % (GNSSsystems[i]))


    fid.write( '\n\n\n\nEND OF OUTPUT FILE')
    fid.close()
    return
