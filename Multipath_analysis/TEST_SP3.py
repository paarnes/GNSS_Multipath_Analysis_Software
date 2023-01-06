from read_SP3Nav import readSP3Nav
filename = 'test1.eph'
sat_pos, epoch_dates, navGNSSsystems, nEpochs, epochInterval,success = readSP3Nav(filename)