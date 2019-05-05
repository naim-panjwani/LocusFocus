function plot_heatmap(genes, tissues, SSPvalues) {
    var data = [{
          z: SSPvalues,
          x: genes,
          y: tissues,
          type: 'heatmap'
    }];
    Plotly.newPlot('heatmap', data);
}