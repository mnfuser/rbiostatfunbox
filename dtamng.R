#     DTAMNG
#
#     Data management
#
#     GNU GPLv3
#
#     Copyright (C) 2023  Marco Manfrini
#
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#     info@manfrinistudio.it

# Source ---------------------------------------------
source("/Users/mmanfrini/Code/rbiostatfunbox/rbiostatfunbox/funbox.R")

# PATH ----

setwd("/Users/mmanfrini/Analisi/maritati/PVL2")
wd<-getwd()
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/preprocess")

# LOAD DATASET ----

dset <- read.csv2(paste0(indir,"/dataset.csv"),
                  stringsAsFactors=T)
View(dset)

# RECODE ----

# NEW VARS ----

# FILTER ----

# WRITE OUT DATASETS ----

## Save
