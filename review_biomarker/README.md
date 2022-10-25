# Biomarker Review Exercise

This example leans on a review experience. Imagine you receive a review request for a paper that claims to have found a biomarker signature for the cognitive abilities of human subjects, comprised of three microRNAs measured in blood. Among many other files, the code supplied to you contains the script for validation of the signature, which you find at `scripts/01_script_ANON_original.R`. The code compares the eigenvalue of the expression of that signature in a cohort of human individuals with either mild cognitive impairment (MCI) or normal cognition (control) and finds a significant difference in eigenvalues between the two groups, which the authors of the manuscript take as confirmation that their signature works, which is why they have promptly patented it for clinical use.

## The exercise

It is your job to ascertain that their claims are, in fact, correct. To do this, please examine the supplied data and code for possible errors or inconsistencies that could lead to a spurious finding. Please first try to make sense of the data and code on your own. After you have done this in a manner that you feel exhausted your capabilities, you can take a look at the second script file, containing hints for some of the possible problems. There are also two scripts containing solutions for two of the problems, but I encourage you to try and identify the problems on your own terms, as running the solution is trivially quick and probably will not lead to much of a learning experience.

## How to start

This subfolder contains an R project file (`review_biomarker.Rproj`). Opening this file with RStudio automatically takes care of the working directory and should allow running the code without any changes. First, clone the parent repository (`git clone <url>`), if you have not already done so, and open the R project file. Then, you can open the script files and run the code. Take care to install all R packages that you don't already have from the corresponding repository (CRAN or Bioconductor), using the base 'install.packages()' function for CRAN and the 'BiocManager' for Bioconductor packages (see https://www.bioconductor.org/install/).
