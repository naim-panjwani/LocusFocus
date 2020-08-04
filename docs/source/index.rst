####################################
Documentation for LocusFocus
####################################

`LocusFocus <https://locusfocus.research.sickkids.ca>`_ is a web application to facilitate the exploration of 
a GWAS signal at a particular locus of the genome and its degree of colocalization
with any SNP-level association data (e.g. expression quantitative trait loci for genes within +/- 1Mbp
in the relevant GTEx tissues selected). 

When paired with GTEx data, the aim is to annotate a GWAS (or region-based association)
to the most probable gene(s) and tissue(s) that may be driving the observed GWAS signal.

In addition, users may upload other datasets to test colocalization with. 
For example, other phenotypic associations (i.e. PheWAS) may be uploaded for assessing pleiotropy, 
or eQTL data from other sources to obtain a formal colocalization test and visualization of the data.

The `Simple Sum method <https://doi.org/10.1371/journal.pgen.1008007>`_
is used for assessing the degree of colocalization of any two given datasets.
When applied to GTEx, LocusFocus presents the degree of colocalization of genes nearby 
the GWAS association for all the tissues selected in an interactive heatmap plot.  

COLOC2 colocalization testing is also available, 
and more colocalization methods may be made available in future version releases.


.. toctree::
   :maxdepth: 3
   :caption: Table of Contents:

   quick_start
   session_retrieval
   examples
   current_versions
   local_installation
   help
   license

