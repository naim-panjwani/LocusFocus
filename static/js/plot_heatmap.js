function plot_heatmap(genes, tissues, SSPvalues) {
    var data = [{
          z: SSPvalues,
          x: genes.map(gene => `<i>${gene}</i>`),
          y: tissues,
          colorscale: 'Portland',
          type: 'heatmap',
          name: '-log10(Simple Sum P-value)',
          hovertemplate: 'Gene: %{x}' +
                         '<br>Tissue: %{y}<br>' +
                         '-log10(SS P-value): %{z}',
          colorbar: {title: '-log10(Simple<br>Sum P-value)'}
    }];
    var layout = {
        // annotations: [],
        margin: {
            r: 50,
            t: 50,
            b: 80,
            l: 300
        }
    }

    // Tried to add the SSPvalue numbers, but does not place correctly (they all go into the middle of the plot)
    // for ( var i = 0; i < tissues.length; i++ ) {
    //     for ( var j = 0; j < genes.length; j++ ) {
    //       var currentValue = SSPvalues[i][j];
    //       if (currentValue != -1) {
    //         var textColor = 'white';
    //       }else{
    //         var textColor = 'black';
    //       }
    //       var result = {
    //         x: genes[j],
    //         y: tissues[i],
    //         text: SSPvalues[i][j],
    //         showarrow: false,
    //         font: {
    //           color: textColor
    //         }
    //       };
    //       layout.annotations.push(result);
    //     }
    // }



    Plotly.newPlot('heatmap', data, layout);
}