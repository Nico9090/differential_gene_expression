tar -xvzf NHC_I2_DIA_IO_L15_Slot1-21_1_4454.d.tar.gz
cd NHC_I2_DIA_IO_L15_Slot1-21_1_4454.d
head analysis.tdf
cd ../
msconvert "NHC_I2_DIA_IO_L15_Slot1-21_1_4454.d" --mzML --outdir "path/to/outDir"
