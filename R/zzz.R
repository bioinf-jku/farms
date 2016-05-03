.onLoad <- function(libname, pkgname) {


	
	require("methods", quietly=TRUE)
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("farms")
	
	
		
		packageStartupMessage(
				" _                             " ,"\n",
				"| |                            " ,"\n",
				"| |  __,   ,_    _  _  _    ,  " ,"\n",
				"|/  /  |  /  |  / |/ |/ |  / \\_" ,"\n",
				"|__/\\_/|_/   |_/  |  |  |_/ \\/ " ,"\n",
				"|\\                             " ,"\n",
				"|/   " ,"\n")	
		
		
		
		
		
	
	
    require(utils)
	require(Biobase, quietly=TRUE) || stop("cannot load farms without Biobase")
	require(affy, quietly=TRUE) || stop("cannot load farms without affy")
    version <- packageDescription("farms",fields="Version")
    packageStartupMessage("Citation: S. Hochreiter et al.,","\n",
						  "A new summarization method for affymetrix probe level data,","\n",
						  "Bioinformatics, 22, 8, 943-949, 2006","\n","\n",
						  "Citation: W. Talloen et al.,","\n",
						  "I/NI-calls for the exclusion of non-informative genes: a highly effective filtering tool for microarray data,","\n",
						  "Bioinformatics, 23, 21, 2897-2902, 2007","\n",
						  "BibTex: enter 'toBibtex(citation(\"farms\"))'","\n\n",
						  "Homepage: http://www.bioinf.jku.at/software/farms/farms.html","\n\n",
						  "FARMS Package Version ", version, "\n")
	
	
	packageStartupMessage("\n","Changes in FARMS:","\n",
						  "For all changes previous to 1.3.0, see the farms vignette.","\n",
						  "Version 1.3.0: Added I/NI-calls for filtering","\n",
						  "               Adjusted Hyperparameters for alternative CDFs,","\n",
						  "               probes set standardized, weighted mean", "\n",
						  "               Works now with R >= 2.8 and Bioconductor 2.3,","\n",
						  "               Changed termination criterion, initialization values,", "\n",
						  "               factors and loadings scaled, added argument robust","\n",
						  "               Update for R-2.11","\n",
						  "               Updated I/NI-Call for Laplace-FARMS version,","\n",
						  "               Maximum likelihood correlation structure given","\n",
						  "               non-negative constraints", "\n",
						  "Version 1.4.0: Default centering changed to median","\n",
						  "Version 1.8.x: Suppression of spurious correlation (Laplace-FARMS)","\n")
	
	upDate.generateExprSet.methods(c(generateExprSet.methods(), "farms"))
	upDate.express.summary.stat.methods(c(express.summary.stat.methods(), "farms"))
	


}
