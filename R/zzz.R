.onLoad <- function(libname, pkgname) {


	
	require("methods", quietly=TRUE)
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("farms")
	
	
	
	packageStartupMessage("#######    #    ######  #     #  ##### " ,"\n",
						  "#         # #   #     # ##   ## #     #" ,"\n",
						  "#        #   #  #     # # # # # #      " ,"\n",
						  "#####   #     # ######  #  #  #  ##### " ,"\n",
						  "#       ####### #   #   #     #       #" ,"\n",
						  "#       #     # #    #  #     # #     #" ,"\n",
						  "#       #     # #     # #     #  ##### " ,"\n")
	
    require(utils)
	require(Biobase, quietly=TRUE) || stop("cannot load farms without Biobase")
	require(affy, quietly=TRUE) || stop("cannot load farms without affy")
    version <- packageDescription("farms",fields="Version")
    packageStartupMessage("Citation: S. Hochreiter et al.,","\n",
						  "A new summarization method for affymetrix probe level data,","\n",
						  "Bioinformatics, 22, 8, 943-949, 2010","\n",
						  "BibTex: enter 'toBibtex(citation(\"farms\"))'","\n\n",
						  "Homepage: http://www.bioinf.jku.at/software/farms/farms.html","\n\n",
						  "FARMS Package Version ", version, "\n")
	
	
	packageStartupMessage("\n","Changes in FARMS:","\n",
						  "For all changes previous to 1.3.0, see the farms vignette.","\n",
						  "Version 1.3.0: Added I/NI-calls for filtering","\n",
						  "Version 1.3.1: Adjusted Hyperparameters for alternative CDFs,","\n",
						  "				probes set standardized, weighted mean", "\n",
						  "Version 1.4.0: Works now with R >= 2.8 and Bioconductor 2.3,","\n",
						  "				Changed termination criterion, initialization values,", "\n",
						  "				factors and loadings scaled, added argument robust","\n",
						  "Version 1.4.1: Update for R-2.11","\n",
						  "Version 1.5.0: Updated I/NI-Call for Laplace-FARMS version,","\n",
						  "				Maximum likelihood correlation structure given","\n",
						  "				non-negative constraints", "\n")
	
	upDate.generateExprSet.methods(c(generateExprSet.methods(), "farms"))
	upDate.express.summary.stat.methods(c(express.summary.stat.methods(), "farms"))
	


}
