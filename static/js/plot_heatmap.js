function plot_heatmap(data, genesdata, SSPvalues) {
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
    // console.log(genenames);
    // console.log(SSPvalues);

    var data = [{
          z: [[1, 20, 30, 50, 1], [20, 1, 60, 80, 30], [30, 60, 1, -10, 20]],
          x: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
          y: ['Morning', 'Afternoon', 'Evening'],
          type: 'heatmap'
    }];
      
    Plotly.newPlot('myDiv', data);
}