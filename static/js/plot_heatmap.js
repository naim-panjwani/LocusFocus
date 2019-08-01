function plot_heatmap(genes, tissues, SSPvalues) {
    var data = [{
          z: SSPvalues,
          x: genes,
          y: tissues,
          type: 'heatmap',
          name: '-log10(Simple Sum P-value)',
          hovertemplate: 'Gene: %{x}' +
                         '<br>Tissue: %{y}<br>' +
                         '-log10(SS P-value): %{z}',
          colorbar: {title: '-log10(Simple<br>Sum P-value)'}
    }];
    var layout = {
        margin: {
            r: 50,
            t: 50,
            b: 80,
            l: 225
        }
    }
    Plotly.newPlot('heatmap', data, layout);
}