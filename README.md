Usage: `process-otu-data.py INPUT-FILE OUTPUT-FILE`

This script processes files containing data on relative abundance of OTUs in one or more samples.

Input should be a CSV file with rows representing different OTUs and columns representing the sampling site(s). The first column should be labeled "OTU" and have the names of each OTU. One of the columns should be labeled "control" and have the relative abundance data for the elution buffer. All other columns can be labeled freely and should have the relative abundance data for each sampling site. Output is another CSV file with columns for each taxonomic rank followed by columns for each sampling site. Each row represents the organism that an OTU is attributed to and each numerical value represents the relative abundance of the OTU converted to a proportion after pruning.

Pruning is done according to relative abundance in the elution buffer. If an OTU has less than PRUNING_VALUE times as many reads in a sample as in the buffer, its presence is attributed to contamination and the number of hits is adjusted to zero. PRUNING_VALUE has been set to 1.5 by default, but may be changed.

Open Tree of Life API code is based on https://github.com/brunoasm/TaxReformer
