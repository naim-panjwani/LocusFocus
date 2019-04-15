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
    ld_colors = [];
    var regionsize = endbp - startbp;
    var log10pvalue_range = d3.max(pvalues.map(p => -Math.log10(p))) - d3.min(pvalues.map(p => -Math.log10(p)));
    var extra_x_range = 0.05 * regionsize;
    var extra_y_range = 0.05 * log10pvalue_range;
    var eqtl_smoothing_window_size = (regionsize/100000) * 15;

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

      while(curr < d3.min([xrange[1],x[x.length-1]]) && indices.length>0) {
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
      }
      else {
        geneRows[currRow].push(geneBin);
      }
    }
    console.log(genesdata);
    console.log(geneRows);
    var gene_area_height = d3.min([0.1 * log10pvalue_range * geneRows.length, 0.5 * log10pvalue_range]);
    var gene_margin = (gene_area_height / geneRows.length) * 0.05;
    var gene_height = (gene_area_height / geneRows.length) - gene_margin;

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
          if(k === 'seq_region_start')  {
            gtex_positions[gtex_tissues[i]].push(+eqtl[k]);
          }
          else if (k === 'minus_log10_p_value') {
            gtex_log10_pvalues[gtex_tissues[i]].push(+eqtl[k]);
          }
          else if (k === 'snp') {
            gtex_snps[gtex_tissues[i]].push(eqtl[k]);
          }
        });
      });
      gtex_line_traces[gtex_tissues[i]] = smoothing(gtex_positions[gtex_tissues[i]], gtex_log10_pvalues[gtex_tissues[i]], 
          [startbp, endbp], eqtl_smoothing_window_size);
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

    var gwas_ymax = d3.max(log10pvalues);
    var gtex_ymax = getYmax(gtex_line_traces);

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
      zeroline: true,
      hovermode: "closest"
    };

    // Plot the genes
    

  Plotly.newPlot('plot', all_traces, layout);
}

