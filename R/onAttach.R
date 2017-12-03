.onAttach <- function(libname, pkgname)
{
    packageStartupMessage("\nThe package 'chngpt' has been loaded. If used for a publication, please cite:\n")
    packageStartupMessage("    Fong, Huang, Gilbert, Permar (2017), BMC Bioinformatics, DOI:10.1186/s12859-017-1863-x")
    packageStartupMessage("    chngpt: threshold regression model estimation and inference\n")
    
#    cat("\n")
#    cat("'chngpt' has been loaded.\n\n")
#    cat("for references type 'citation()' and 'citation('chngpt')'.\n\n")
}
