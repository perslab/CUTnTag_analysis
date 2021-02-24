
R -e "rmarkdown::render('03_demultiplex.Rmd',output_file='03_demultiplex.html')"
bash -x 04_filteredCells_makeBAMmakeBWfindPeaks.sh &> logfile_04.log
R -e "rmarkdown::render('05_Bin5000_Processing.Rmd',output_file='05_Bin5000_Processing.html')"
bash -x 06_BinClusters_makeBAMmakeBWfindPeaks.sh &> logfile_06.log
R -e "rmarkdown::render('07_peakMatrix_processing.Rmd',output_file='07_peakMatrix_processing.html')"
bash -x 08_peakClusters_makeBAMmakeBWfindPeaks.sh &> logfile_08.log
