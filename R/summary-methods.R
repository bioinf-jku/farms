##    Copyright (C) 2006 Djork-Arne Clevert (okko@clevert.de),
##       				 Sepp Hochreiter (hochreit@bioinf.jku.at),
##                       Klaus Obermayer (oby@cs.tu-berlin.de)
##    Berlin University of Technology,
##    Institute for Software Engineering and Theoretical Computer Science 
##    The software is maintained and developed by Djork-Arné Clevert. 
##    We offer a first implementation of the new 
##    ``Factor Analysis for Robust Microarray Summarization'' (FARMS) algorithm.
##    This program is free software; you can redistribute it and/or modify it under 
##    the terms of the GNU General Public License as published by the Free Software 
##    Foundation; either version 2 of the License, or (at your option) any later version. 
##    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
##    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
##    See the GNU General Public License for more details.
##    If you use this library, please cite:
##
##    @article{SeppHochreiter02102006,
##		author = {Hochreiter, Sepp and Clevert, Djork-Arne and Obermayer, Klaus},
##		title = {{A new summarization method for Affymetrix probe level data}},
##		journal = {Bioinformatics},
##		volume = {},
##		number = {},
##		pages = {btl033},
##		doi = {10.1093/bioinformatics/btl033},
##		year = {2006},
##		URL = {http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl033v1},
##		eprint = {http://bioinformatics.oxfordjournals.org/cgi/reprint/btl033v1.pdf}
##		}

##	@article{INI-Calls:07,
##		author = {Willem Talloen and Djork-Arne Clevert and Sepp Hochreiter and Dhammika Amaratunga and Luc Bijnens and Stefan Kass and Hinrich W.H. Ghlmann},
##		title = {/NI-calls for the exclusion of non-informative genes: a highly effective filtering tool for microarray data},
##		journal = {Bioinformatics},
##		volume = {},
##		number = {},
##		pages = {btm478},
##		doi = {doi:10.1093/bioinformatics/btm478},
##		year = {2007},
##		URL = {http://bioinformatics.oxfordjournals.org/cgi/content/short/btm478v1},
##		eprint = {http://bioinformatics.oxfordjournals.org/cgi/reprint/btm478v1}
##		}





setMethod("summary","INI_Calls",function(object,...){
   cat("Summary \n")
   cat("Informative probe sets      : ",round(100 * length(object@I_Calls)/(length(object@I_Calls)+length(object@NI_Calls)),digits=2),"% \n",sep="")
   cat("Non-Informative probe sets  : ",round(100 * length(object@NI_Calls)/(length(object@I_Calls)+length(object@NI_Calls)),digits=2),"% \n",sep="")	# number of informative genes
})

