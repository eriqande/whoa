# whoa 0.0.1

* First submission to CRAN
* Added a `NEWS.md` file to track changes to the package.



# whoa 0.0.2 2018-11-09

* `read_whoa`: new function that reads VCF file for whoa using SeqArray
* `run_whoa`: new wrapper function that runs all the functions necessary to 
generate the plots and the heterozygote miscall rate
* modified existing codes to use SeqArray (vcfR is still there).
* vcfR and SeqArray are in the Suggests field and checks that packages are
installed depending on class of objet used.
