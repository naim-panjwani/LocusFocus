function plot_heatmap(data, genesdata, pvaluedata) {
    var positions = data.positions;
    var pvalues = data.pvalues;
    var ld_values = data.ld_values;
    var chrom = data.chrom;
    if(chrom === 23) chrom = "X";
    var startbp = data.startbp;
    var endbp = data.endbp;
    var snps = data.snps;
    var lead_snp = data.lead_snp;
    var gtex_tissues = data.gtex_tissues;

    var genenames = genesdata.map(gene => gene.name);
    console.log(genenames);
}