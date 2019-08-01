function plot_heatmap(genes, tissues, SSPvalues) {
    var data = [{
          z: SSPvalues,
          x: genes,
          y: tissues,
          type: 'heatmap',
          hovertemplate: 'Gene: %{x}' +
                         '<br>Tissue: %{y}<br>' +
                         '-log10(SS P-value): %{z}'
    }];
    Plotly.newPlot('heatmap', data);
}