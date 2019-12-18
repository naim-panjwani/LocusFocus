function copy(aObject) {
    if (!aObject) {
        return aObject;
    }

    let v;
    let bObject = Array.isArray(aObject) ? [] : {};
    for (const k in aObject) {
        v = aObject[k];
        bObject[k] = (typeof v === "object") ? copy(v) : v;
    }

    return bObject;
}

function plot_heatmap(genes, tissues, SSPvalues, image_width=1080, image_height=1080) {
    // remember that these are -log10P-values
    // want the different negative p-value statuses to have white/black/grey colors
    // can't really do it that way b/c of the different cases that may occur
    // prior checks ensure that we have at least one positive -log10 SS Pvalue
    var pmax=0;
    var newSSPvalues = copy(SSPvalues);
    for(i=0; i<tissues.length; i++) {
        for(j=0; j<genes.length; j++) {
            newSSPvalues[i][j] = parseFloat(newSSPvalues[i][j]);
            if(newSSPvalues[i][j] === -2) {
                newSSPvalues[i][j] = -1;
            } else if(newSSPvalues[i][j] === -3) {
                newSSPvalues[i] = -1;
            } else if(newSSPvalues[i][j] > pmax) {
                pmax = newSSPvalues[i][j];
            }
        }
    }
    var colorscale_exception_percentage = 0.05;
    var new_minp = -1 * (colorscale_exception_percentage / (1-colorscale_exception_percentage)) * pmax;
    for(i=0; i<tissues.length; i++) {
        for(j=0; j<genes.length; j++) {
            if(newSSPvalues[i][j] === -1) {
                newSSPvalues[i][j] = new_minp;
            }
        }
    }
    var data = [{
          z: newSSPvalues,
          x: genes.map(gene => `<i>${gene}</i>`),
          y: tissues,
          //colorscale: 'Portland',
          type: 'heatmap',
          name: '-log10(Simple Sum P-value)',
          hovertemplate: 'Gene: %{x}' +
                         '<br>Tissue: %{y}<br>' +
                         `-log10(SS P-value): %{z}`,
          colorbar: {title: '-log10(Simple<br>Sum P-value)', dtick0: 0, dtick: 1, autotick: false},
          colorscale: [
            // Values between 0-15% of the min and max of z; ie. our negative stasuses

            //   [0, 'rgb(0,0,0)'], // originally -3; black
            //   [0.05, 'rgb(0,0,0)'], // first 5% is black
              
            //   [0.05, 'rgb(255,255,255)'], // originally -2; white
            //   [0.10, 'rgb(255,255,255)'], // between 5-10%

              [0, 'rgb(105,105,105)'], // originally -1; gray
              [0.05, 'rgb(105,105,105)'], 

              // Next trying to follow LD colors

              // dark blue to bright red HSV gradient in 5 steps:
              [0.05, 'rgb(0, 0, 128)'], // dark blue
              [0.24, 'rgb(0, 0, 128)'], // dark blue
              
              [0.24, 'rgb(0, 147, 142)'],
              [0.43, 'rgb(0, 147, 142)'],
              
              [0.43, 'rgb(10, 166, 0)'],
              [0.62, 'rgb(10, 166, 0)'],
              
              [0.62, 'rgb(185, 167, 0)'],
              [0.81, 'rgb(185, 167, 0)'],
              
              [0.81, 'rgb(204, 0, 24)'],
              [1.0, 'rgb(204, 0, 24)'] // using darker red instead of rgb(255,0,0) bright red

              // // Values between 20-30% of the min and max of z
              // // have color rgb(40, 40, 40)

              // [0.2, 'rgb(40, 40, 40)'],
              // [0.3, 'rgb(40, 40, 40)'],

              // [0.3, 'rgb(60, 60, 60)'],
              // [0.4, 'rgb(60, 60, 60)'],

              // [0.4, 'rgb(80, 80, 80)'],
              // [0.5, 'rgb(80, 80, 80)'],

              // [0.5, 'rgb(100, 100, 100)'],
              // [0.6, 'rgb(100, 100, 100)'],

              // [0.6, 'rgb(120, 120, 120)'],
              // [0.7, 'rgb(120, 120, 120)'],

              // [0.7, 'rgb(140, 140, 140)'],
              // [0.8, 'rgb(140, 140, 140)'],

              // [0.8, 'rgb(160, 160, 160)'],
              // [0.9, 'rgb(160, 160, 160)'],

              // [0.9, 'rgb(180, 180, 180)'],
              // [1.0, 'rgb(180, 180, 180)']
          ]
    }];
    var layout = {
        // annotations: [],
        margin: {
            r: 50,
            t: 50,
            b: 125,
            l: 300
        },
        autosize: false,
        width: image_width,
        height: image_height
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


