.. _future:

######################################
Future Directions
######################################

LocusFocus is `open source <https://github.com/naim-panjwani/LocusFocus>`_ and continuously undergoing development and improvements.  

Below is a list of improvements and features currently under development:

- Enable uploading of compressed files (then convert to bgzip/tabix)
- Enable immediate deletion of the session data
- Enable option to match variant names by just chrom:pos information (currently should only be either all rsid's or all variant_id format)
- Generate updated 1KG Phase 3 binary PLINK files for hg19 LD calculations (currently only have 2012 1KG file)
- Enable user to set the window sizes for eQTL lines; and make calculation clearer
- Make a table of GWAS and eQTL merged results
- Plot of beta correlations
- Make P-P plot


A list of known bugs currently being addressed:

- chrX plots do not show correctly
- SS coordinate input field checking not working - should also check if it's a subregion of full coordinate/locus region
- Remove related individuals prior to calculating LD from the 1000 Genomes
- Include option for NFE LD (Non-Finnish European subset) for hg19 1000 Genomes
- File too large error handler not working
- Fix COLOC2 chromosome X issue
- Add requirement for having REF and ALT columns for secondary datasets for better SNP matching
