function plot_gwas(data, genesdata) {
  
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
    var SS_start = +data.SS_region[0];
    var SS_end = +data.SS_region[1];

    // console.log(data);
    // console.log(genesdata);

    var log10pvalues = pvalues.map(p => -Math.log10(p))
    var no_ld_info_snps = [], no_ld_info_snps_color = "#7f7f7f"; // grey
    var ld_lt_20_group = [], ld_lt_20_group_color = "#1f77b4"; // very low ld points (dark blue)
    var ld_20_group = [], ld_20_group_color = "#17becf"; // low ld points (light blue)
    var ld_40_group = [], ld_40_group_color = "#bcbd22"; // green
    var ld_60_group = [], ld_60_group_color = "#ff7f0e"; // orange
    var ld_80_group = [], ld_80_group_color = "#d62728"; // red
    var lead_snp_index = 0, lead_snp_color = "#9467bd"; // purple
    var markersize = 10;
    var lead_markersize = 10*1.5;
    var ld_colors = [];
    var regionsize = endbp - startbp;
    var log10pvalue_range = d3.max(pvalues.map(p => -Math.log10(p))) - d3.min(pvalues.map(p => -Math.log10(p)));
    var extra_x_range = 0.05 * regionsize;
    var extra_y_range = 0.05 * log10pvalue_range;
    var eqtl_smoothing_window_size = (regionsize/100000) * 15;
    var percent_occupied_by_one_char = 0.011;
    var font_height = 0.5;
    // console.log(eqtl_smoothing_window_size)

    // Helper functions:
    function smoothing(x,y,xrange,window_size) {
      window_partition = window_size;
      windowing = (xrange[1]-xrange[0])/window_partition;
      curr = xrange[0];
      smooth_curve_x = [];
      smooth_curve_y = [];
      indices = [];
      x.map((v,i) => {
        if(v>=curr && v<=(curr+windowing)) indices.push(i);
      });

      while(curr < d3.min([xrange[1],x[x.length-1]])) {
        if(indices.length > 0) {
          yAtIndices = indices.map(i => y[i]);
          xAtIndices = indices.map(i => x[i]);
          ymaxAtIndices = d3.max(yAtIndices);
          if(ymaxAtIndices === -1 || ymaxAtIndices < 0) ymaxAtIndices=0;
          desiredXindex = [];
          indices.map(i => {
            if(y[i] === ymaxAtIndices) desiredXindex.push(i);
          });
          smooth_curve_x.push(x[desiredXindex[0]]);
          smooth_curve_y.push(ymaxAtIndices);
          curr = curr + windowing + 1;
          indices = [];
          x.map((v,i) => {
            if(v>=curr && v<=(curr+windowing)) indices.push(i);
          });
        }
        else {
          curr = curr + windowing + 1;
          indices = [];
          x.map((v,i) => {
            if(v>=curr && v<=(curr+windowing)) indices.push(i);
          });
        }
      }
      return [smooth_curve_x, smooth_curve_y];
    }

    function getYmax(gtex_traces) {
      ymax = 1;
      for(var i=0; i<gtex_tissues.length; i++) {
        currymax = d3.max(gtex_traces[gtex_tissues[i]][1]);
        if(currymax > ymax) ymax = currymax;
      }
      return ymax;
    }
    

    // Add the row number each gene will be plotted into:
    genesdata.sort(function (a, b) {
      return a.txStart - b.txStart;
    });

    function overlap(bin1, bin2) {
      return (bin2[0]>=bin1[0] && bin2[0]<=bin1[1]) || (bin2[1]>=bin1[0] && bin2[1]<=bin1[1]);
    }

    function checkSpace(geneRows, rowNum, bin) {
      occupied = false;
      for(var i=0; i<geneRows[rowNum].length; i++) {
        if(overlap(geneRows[rowNum][i], bin)) {
          occupied = true;
          return occupied;
        }
      }
      return occupied;
    }

    var firstBin = [genesdata[0].txStart, genesdata[0].txEnd];
    var geneRows = [];
    geneRows[0] = [firstBin];
    genesdata[0]['geneRow'] = 1;

    for(var i = 1; i < genesdata.length; i++) {
      currRow = 0;
      geneBin = [genesdata[i].txStart, genesdata[i].txEnd];
      occupied = checkSpace(geneRows, currRow, geneBin);
      while(occupied && currRow < geneRows.length) {
        currRow++;
        if(currRow < geneRows.length) occupied = checkSpace(geneRows, currRow, geneBin);
      }
      if(currRow === geneRows.length) {
        geneRows[currRow] = [geneBin];
        genesdata[i]['geneRow'] = currRow+1;
      }
      else {
        geneRows[currRow].push(geneBin);
        genesdata[i]['geneRow'] = currRow+1;
      }
    }

    // console.log(genesdata);
    // console.log(geneRows);
    var gene_area_percentage = 0.15;
    var max_num_gene_rows = 7;
    var gene_area_height = d3.min([gene_area_percentage * log10pvalue_range * geneRows.length, gene_area_percentage * max_num_gene_rows * log10pvalue_range]);
    var row_height = gene_area_height / geneRows.length;
    var text_height = row_height * 0.15
    var gene_margin = row_height * 0.15;
    var exon_height = row_height - (2 * (gene_margin + text_height));
    var intron_height = exon_height * 0.4;

    // console.log(row_height);
    // console.log(text_height*2+gene_margin*2+exon_height);

    var rectangle_shapes = [];
    var annotations_x = [];
    var annotations_y = [];
    var annotations_text = [];
    for(var i=0; i < genesdata.length; i++) {
      // build intron rectangle shapes for each gene:
      var rectangle_shape = {
        type: 'rect',
        xref: 'x',
        yref: 'y',
        x0: genesdata[i]['txStart'],
        y0: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + ((exon_height - intron_height)/2),
        x1: genesdata[i]['txEnd'],
        y1: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + ((exon_height - intron_height)/2) + intron_height,
        line: {
          color: 'rgb(55, 128, 191)',
          width: 1
        },
        fillcolor: 'rgba(55, 128, 191, 1)'
      }
      rectangle_shapes.push(rectangle_shape);
      annotations_x.push(genesdata[i]['txStart']);
      annotations_x.push((genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2);
      annotations_x.push(genesdata[i]['txEnd']);
      var y = -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + ((exon_height - intron_height)/2) + intron_height/2;
      annotations_y.push(y);
      annotations_y.push(y);
      annotations_y.push(y);
      annotations_text.push(genesdata[i]['name']);
      annotations_text.push(genesdata[i]['name']);
      annotations_text.push(genesdata[i]['name']);
      for(var j=0; j < genesdata[i]['exonStarts'].length; j++) {
        // build exon rectangle shapes for current gene
        var rectangle_shape = {
          type: 'rect',
          xref: 'x',
          yref: 'y',
          x0: genesdata[i]['exonStarts'][j],
          y0: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin,
          x1: genesdata[i]['exonEnds'][j],
          y1: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + exon_height,
          line: {
            color: 'rgb(55, 128, 191)',
            width: 1
          },
          fillcolor: 'rgba(55, 128, 191, 1)'
        }
        rectangle_shapes.push(rectangle_shape);
      }
    }

    // Smooth out each GTEx tissue's association results for plotting as lines:
    gtex_line_traces = {};
    gtex_positions = {};
    gtex_log10_pvalues = {};
    gtex_snps = {};
    for(var i = 0; i < gtex_tissues.length; i++) {
      gtex_positions[gtex_tissues[i]] = [];
      gtex_log10_pvalues[gtex_tissues[i]] = [];
      gtex_snps[gtex_tissues[i]] = [];
      data[gtex_tissues[i]].forEach(eqtl => {
        Object.keys(eqtl).forEach(k => {
          if(k === 'variant_pos')  {
            gtex_positions[gtex_tissues[i]].push(+eqtl[k]);
          }
          else if (k === 'pval') {
            gtex_log10_pvalues[gtex_tissues[i]].push(-Math.log10(+eqtl[k]));
          }
          else if (k === 'rs_id') {
            gtex_snps[gtex_tissues[i]].push(eqtl[k]);
          }
        });
      });
      gtex_line_traces[gtex_tissues[i]] = smoothing(gtex_positions[gtex_tissues[i]], gtex_log10_pvalues[gtex_tissues[i]], 
          [startbp, endbp], eqtl_smoothing_window_size);
      // console.log(gtex_line_traces);
    }    

    // Assign each SNP to an LD group:
    for(i=0; i<ld_values.length; i++) {
      if (snps[i] === lead_snp) {
        lead_snp_index = i;
        ld_colors[i] = lead_snp_color;
      }
      else if (ld_values[i] === -1 || ld_values[i] < 0) { 
       no_ld_info_snps.push(i);
       ld_colors[i] = no_ld_info_snps_color;
      }
      else if (Math.abs(ld_values[i]) < 0.2) {
        ld_lt_20_group.push(i);
        ld_colors[i] = ld_lt_20_group_color;
      }
      else if (Math.abs(ld_values[i]) >= 0.2 && Math.abs(ld_values[i]) < 0.4) {
        ld_20_group.push(i);
        ld_colors[i] = ld_20_group_color;
      }
      else if (Math.abs(ld_values[i]) >= 0.4 && Math.abs(ld_values[i]) < 0.6) {
        ld_40_group.push(i);
        ld_colors[i] = ld_40_group_color;
      }
      else if (Math.abs(ld_values[i]) >= 0.6 && Math.abs(ld_values[i]) < 0.8) {
        ld_60_group.push(i);
        ld_colors[i] = ld_60_group_color;
      }
      else if (Math.abs(ld_values[i]) >= 0.8) {
        ld_80_group.push(i);
        ld_colors[i] = ld_80_group_color;
      }
      else if (snps[i] == lead_snp) {
        lead_snp_index = i;
        ld_colors[i] = lead_snp_color;
      }
      else {
        no_ld_info_snps.push(i);
        ld_colors[i] = no_ld_info_snps_color;
      }
    }
    
    // plot the 7 LD groups
    var no_ld_trace = {
      x: no_ld_info_snps.map(i => positions[i]),
      y: no_ld_info_snps.map(i => log10pvalues[i]),
      name: 'No LD Info',
      mode: 'markers',
      type: 'scatter',
      text: no_ld_info_snps.map(i => snps[i]),
      marker: {
        size: markersize,
        color: no_ld_info_snps_color
      }
    };

    var ld_lt_20_trace = {
      x: ld_lt_20_group.map(i => positions[i]),
      y: ld_lt_20_group.map(i => log10pvalues[i]),
      name: '< 0.2',
      mode: 'markers',
      type: 'scatter',
      text: ld_lt_20_group.map(i => snps[i]),
      marker: {
        size: markersize,
        color: ld_lt_20_group_color
      }
    };

    var ld_20_trace = {
      x: ld_20_group.map(i => positions[i]),
      y: ld_20_group.map(i => log10pvalues[i]),
      name: '0.2',
      mode: 'markers',
      type: 'scatter',
      text: ld_20_group.map(i => snps[i]),
      marker: {
        size: markersize,
        color: ld_20_group_color
      }
    };
  
    var ld_40_trace = {
      x: ld_40_group.map(i => positions[i]),
      y: ld_40_group.map(i => log10pvalues[i]),
      name: '0.4',
      mode: 'markers',
      type: 'scatter',
      text: ld_40_group.map(i => snps[i]),
      marker: {
        size: markersize,
        color: ld_40_group_color
      }
    };

    var ld_60_trace = {
      x: ld_60_group.map(i => positions[i]),
      y: ld_60_group.map(i => log10pvalues[i]),
      name: '0.6',
      mode: 'markers',
      type: 'scatter',
      text: ld_60_group.map(i => snps[i]),
      marker: {
        size: markersize,
        color: ld_60_group_color
      }
    };

    var ld_80_trace = {
      x: ld_80_group.map(i => positions[i]),
      y: ld_80_group.map(i => log10pvalues[i]),
      name: '0.8',
      mode: 'markers',
      type: 'scatter',
      text: ld_80_group.map(i => snps[i]),
      marker: {
        size: markersize,
        color: ld_80_group_color
      }
    };

    var lead_snp_trace = {
      x: [positions[lead_snp_index]],
      y: [log10pvalues[lead_snp_index]],
      name: 'Lead SNP',
      mode: 'markers',
      type: 'scatter',
      text: lead_snp,
      marker: {
        size: lead_markersize,
        color: lead_snp_color
      },
      yaxis: 'y1'
    };

    all_traces = [ no_ld_trace, ld_lt_20_trace, ld_20_trace, ld_40_trace, ld_60_trace, ld_80_trace, lead_snp_trace ];


    // Plot the GTEx lines (gtex_line_traces):
    for(var i=0; i < gtex_tissues.length; i++) {
      var gtex_tissue_trace = {
        x: gtex_line_traces[gtex_tissues[i]][0],
        y: gtex_line_traces[gtex_tissues[i]][1],
        name: gtex_tissues[i],
        mode: 'lines',
        xaxis: 'x1',
        yaxis: 'y2'
      };
      all_traces.push(gtex_tissue_trace);
    }

    for(var i=0; i < gtex_tissues.length; i++) {
      var gtex_tissue_trace = {
        x: gtex_positions[gtex_tissues[i]],
        y: gtex_log10_pvalues[gtex_tissues[i]],
        name: gtex_tissues[i],
        mode: 'markers',
        type: 'scatter',
        text: gtex_snps[gtex_tissues[i]],
        marker: {
          size: markersize,
          opacity: 0.3
        },
        xaxis: 'x1',
        yaxis: 'y2',
        visible: 'legendonly'
      };
      all_traces.push(gtex_tissue_trace);
    }

    var genenames_trace = {
      x: annotations_x,
      y: annotations_y,
      text: annotations_text,
      type: 'scatter',
      mode: 'markers',
      marker: {
        opacity: 0
      },
      yaxis: 'y1',
      showlegend: false,
      name: 'Gene name'
    }
    all_traces.push(genenames_trace);

    var gwas_ymax = d3.max(log10pvalues);
    var gtex_ymax = getYmax(gtex_line_traces);
    
    // Find a place for gene names text

    function point_overlap(point, rect) {
      xcheck = (point[0]>=rect[0] && point[0]<=rect[2]);
      ycheck = (point[1]>=rect[1] && point[1]<=rect[3]);
      return(xcheck && ycheck);
    }

    function rect_overlap(rect1, rect2) {
      var overlap = false;
      if( point_overlap([rect2[0],rect2[1]], rect1) || 
          point_overlap([rect2[3], rect2[2]], rect1) || 
          point_overlap([rect2[0], rect2[3]], rect1) || 
          point_overlap([rect2[2], rect2[3]], rect1)) {
        overlap = true;
      }
      return(overlap);
    }

    function curr_rect_overlaps(rect_bin, rect_bins) {
      overlap = false;
      for(var i=0; i<rect_bins.length; i++) {
        var curr_rect_bin = rect_bins[i];
        if(rect_overlap(rect_bin, curr_rect_bin)) {
          return(true);
        }
      }
      return(overlap);
    }

    var full_y_range = gwas_ymax+extra_y_range+gene_area_height;
    i = 0;
    midx = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
    midy = -(genesdata[i]['geneRow'] * row_height) + row_height/2;
    thegenename = genesdata[i]['name'];
    var xrefloc = (midx - (startbp - extra_x_range)) / (regionsize + 2 * extra_x_range);
    var x0 = d3.max([xrefloc - thegenename.length/2 * percent_occupied_by_one_char, 0]);
    var x1 = xrefloc + thegenename.length/2 * percent_occupied_by_one_char;
    var tempx2 = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
    if (x0<0 && x1>0) {
      x1=x1 + (-1*x0); // move forward by x0
      x0=0; // move up to 0
      tempx2 = (x0/2 + x1/2) * (regionsize + extra_x_range*2) + startbp - extra_x_range;
    }
    if(x0<1 && x1>1) {
      x0=1 - (x1-1); // move back by x1-1
      x1=1;
      tempx2 = (x0/2 + x1/2) * (regionsize + extra_x_range*2) + startbp - extra_x_range;
    }
    var y0 = (gene_area_height - (-1*midy) - font_height) / full_y_range;
    var y1 = (gene_area_height - (-1*midy)) / full_y_range;
    var first_rect_bin = [x0,y0,x1,y1];
    //var roughly_half_row_height = text_height + gene_margin + ((exon_height - intron_height)/2) + intron_height/2;

    var rect_bins = [];
    rect_bins = [first_rect_bin];
    var annotations_x2 = [];
    var annotations_y2 = [];
    var annotations_text2 = [];
    var locations = [];
    annotations_x2.push(tempx2);
    annotations_y2.push(-(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin);
    annotations_text2.push(`<i>${genesdata[i]['name']}</i>`);
    locations.push('bottom');

    // var temp_data = [{
    //   'genename': annotations_text2[i],
    //   'x': annotations_x2[i],
    //   'y': annotations_y2[i],
    //   'rect_bin': first_rect_bin,
    //   'all rect_bins': rect_bins,
    //   'location': locations[i]
    // }];

    for(var i=1; i < genesdata.length; i++) {
      midx = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
      midy = -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin;
      thegenename = genesdata[i]['name'];
      var xrefloc = (midx - (startbp - extra_x_range)) / (regionsize + 2 * extra_x_range);
      var x0 = d3.max([xrefloc - annotations_text[i].length/2 * percent_occupied_by_one_char, 0]);
      var x1 = xrefloc + thegenename.length/2 * percent_occupied_by_one_char;
      var tempx2 = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
      if (x0<0 && x1>0) {
        x1=x1 + (-1*x0); // move forward by x0
        x0=0; // move up to 0
        tempx2 = (x0/2 + x1/2) * (regionsize + extra_x_range*2) + startbp - extra_x_range;
      }
      if(x0<1 && x1>1) {
        x0=1 - (x1-1); // move back by x1-1
        x1=1;
        tempx2 = (x0/2 + x1/2) * (regionsize + extra_x_range*2) + startbp - extra_x_range;
      }
      //console.log(font_height);
      //console.log(text_height);
      var y0 = (gene_area_height - (-1*midy) - font_height) / full_y_range;
      var y1 = (gene_area_height - (-1*midy)) / full_y_range;
      var curr_rect_bin = [x0,y0,x1,y1];
      if(curr_rect_overlaps(curr_rect_bin, rect_bins)) {
        // try the top area of the gene
        y0 = (gene_area_height - (-1*midy) - font_height + exon_height) / full_y_range;
        y1 = (gene_area_height - (-1*midy) + exon_height) / full_y_range;
        curr_rect_bin = [x0,y0,x1,y1];
        if(curr_rect_overlaps(curr_rect_bin, rect_bins)) {
          rect_bins.push([-1,-1,-1,-1]); // don't output the genename text then
          annotations_x2.push(-1);
          annotations_y2.push(-1);
          annotations_text2.push(`<i>${genesdata[i]['name']}</i>`);
          locations.push('hidden');
        } else {
          rect_bins.push(curr_rect_bin); // put gene name text at the top of the gene
          annotations_x2.push(tempx2);
          annotations_y2.push(-(genesdata[i]['geneRow'] * row_height) + row_height);
          annotations_text2.push(`<i>${genesdata[i]['name']}</i>`);
          locations.push('top');
        }
      } else {
        // default to putting the gene name text at the bottom of the gene
        rect_bins.push(curr_rect_bin);
        annotations_x2.push(tempx2);
        annotations_y2.push(-(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin);
        annotations_text2.push(`<i>${genesdata[i]['name']}</i>`);
        locations.push('bottom');
      }
    //   temp_data.push({
    //     'genename': annotations_text2[i],
    //     'x': annotations_x2[i],
    //     'y': annotations_y2[i],
    //     'rect_bin': curr_rect_bin,
    //     'all rect_bins': rect_bins,
    //     'location': locations[i]
    //   });
    }
    // console.log(temp_data);

    
    // SOME TEMPORARY TRIAL CODE FOR DETERMINING TEXT BINNING
    /**
    var i = 13;
    var xrefloc = (((genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2) - (startbp - extra_x_range)) / (regionsize + 2 * extra_x_range);
    var leftside = d3.max([xrefloc - annotations_text2[i].length/2 * percent_occupied_by_one_char, 0]);
    var rightside = xrefloc + annotations_text2[i].length/2 * percent_occupied_by_one_char;
    if(rightside>1) {
      rightside=1
    } else if (rightside<0) {
      rightside=0;
    };
    if(leftside>1) {
      leftside=1
    } else if (leftside<0) {
      leftside=0;
    };
  
    var trial_trace = {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: leftside,
      y0: gene_area_height/full_y_range - (-1*annotations_y2[i])/full_y_range - font_height/full_y_range,
      x1: rightside,
      y1: gene_area_height/full_y_range - (-1*annotations_y2[i])/full_y_range,
      fillcolor: 'red',
      opacity: 0.3,
      line: {width: 1}
    }
    rectangle_shapes.push(trial_trace);
    */

    var genenames_trace2 = {
      x: annotations_x2,
      y: annotations_y2,
      text: annotations_text2,
      type: 'scatter',
      mode: 'markers+text',
      marker: {
        opacity: 0
      },
      yaxis: 'y1',
      showlegend: false,
      name: 'Gene name',
      textposition: 'bottom'
    }
    all_traces.push(genenames_trace2);
    
    // Shade the Simple Sum Region
    var SS_shade_shape = {
      type: 'rect',
      xref: 'x',
      yref: 'y',
      x0: SS_start,
      y0: 0,
      x1: SS_end,
      y1: gwas_ymax,
      fillcolor: "#d3d3d3",
      opacity: 0.3,
      line: { width: 0 }
    }
    rectangle_shapes.push(SS_shade_shape);


    var layout = {
      xaxis: {
        range: [startbp - extra_x_range, endbp + extra_x_range],
        zeroline: false,
        title: { text: `Chromosome ${chrom} (hg19)` }
      },
      yaxis: {
        range: [0 - gene_area_height, gwas_ymax + extra_y_range],
        title: { text: `GWAS -log10(p-value)` }
      },
      yaxis2: {
        range: [0 - gene_area_height * (gtex_ymax/gwas_ymax), gtex_ymax + extra_y_range * (gtex_ymax/gwas_ymax)],
        overlaying: 'y',
        anchor: 'x',
        side: 'right',
        showgrid: false,
        title: 'GTEx eQTL -log10(p-value)'
      },
      height: 700,
      width: 960,
      showlegend: true,
      legend: {
        x: 1.1,
        y: 1
      },
      zeroline: true,
      hovermode: "closest",
      shapes: rectangle_shapes
    };

  Plotly.newPlot('plot', all_traces, layout);
}

