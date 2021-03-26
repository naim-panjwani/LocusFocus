var theTable = d3.select("#variants-table");
var theNTable = d3.select("#numSSsnpsUsed-table");
var secondaryTable = d3.select("#secondary-table")

function buildTable(genes, tissues, SSPvalues, transpose = false) {
    
    var thead = d3.select('#variants-table').select("thead");
    var tbody = d3.select("#variants-table").select("tbody");

    // Clear table:
    if ( $.fn.dataTable.isDataTable( '#variants-table' ) ) {
        var mytable = $('#variants-table').DataTable();
        mytable.destroy();
    }

    thead.text("");
    tbody.text("");
    // thead.innerHTML = "";
    // tbody.innerHTML = "";

    // Build column headers:
    if(transpose) {
        thead.append('tr').append('th').attr('class','th-sm').text('Tissue')
        for(i=0; i<genes.length; i++) {
            theTable.select('thead').select('tr')
                .append('th')
                .attr('class', 'th-sm')
                .text(genes[i]);
        }
    }
    else {
        thead.append('tr').append('th').attr('class','th-sm').text('Gene')
        for(i=0; i<tissues.length; i++) {
            theTable.select('thead').select('tr')
                .append('th')
                .attr('class', 'th-sm')
                .text(tissues[i]);
        }
    }
    

    // Add table body:
    if(transpose) {
        for(i=0; i<tissues.length; i++) { // for each tissue
            var row = tbody.append('tr');
            row.append('td').text(tissues[i]);
            for(j=0; j<genes.length; j++) { // for each column
                row.append('td').text(SSPvalues[i][j]);
            }
        }
    }
    else {
        for(i=0; i<genes.length; i++) { // for each tissue
            var row = tbody.append('tr');
            row.append('td').text(genes[i]);
            for(j=0; j<tissues.length; j++) { // for each column
                row.append('td').text(SSPvalues[j][i]);
            }
        }
    }
    
    // Add DataTables functionality:
    varTable = $(document).ready(function () {
        var thedatatable = $('#variants-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'SS_pvalues'
                },
                {
                    extend: 'excel',
                    filename: 'SS_pvalues',
                    messageTop: 'Simple Sum P-values'
                }
            ]
        });
    });
    
}





function buildNTable(genes, tissues, num_SS_snps_used, transpose = false) {
    
    var thead = d3.select('#numSSsnpsUsed-table').select("thead");
    var tbody = d3.select("#numSSsnpsUsed-table").select("tbody");

    // Clear table:
    if ( $.fn.dataTable.isDataTable( '#numSSsnpsUsed-table' ) ) {
        var mytable = $('#numSSsnpsUsed-table').DataTable();
        mytable.destroy();
    }

    thead.text("");
    tbody.text("");
    // thead.innerHTML = "";
    // tbody.innerHTML = "";

    // Build column headers:
    if(transpose) {
        thead.append('tr').append('th').attr('class','th-sm').text('Tissue')
        for(i=0; i<genes.length; i++) {
            theNTable.select('thead').select('tr')
                .append('th')
                .attr('class', 'th-sm')
                .text(genes[i]);
        }
    }
    else {
        thead.append('tr').append('th').attr('class','th-sm').text('Gene')
        for(i=0; i<tissues.length; i++) {
            theNTable.select('thead').select('tr')
                .append('th')
                .attr('class', 'th-sm')
                .text(tissues[i]);
        }
    }
    

    // Add table body:
    if(transpose) {
        for(i=0; i<tissues.length; i++) { // for each tissue
            var row = tbody.append('tr');
            row.append('td').text(tissues[i]);
            for(j=0; j<genes.length; j++) { // for each column
                row.append('td').text(num_SS_snps_used[i][j]);
            }
        }
    }
    else {
        for(i=0; i<genes.length; i++) { // for each tissue
            var row = tbody.append('tr');
            row.append('td').text(genes[i]);
            for(j=0; j<tissues.length; j++) { // for each column
                row.append('td').text(num_SS_snps_used[j][i]);
            }
        }
    }
    
    // Add DataTables functionality:
    varTable2 = $(document).ready(function () {
        var thedatatable2 = $('#numSSsnpsUsed-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'numSSsnpsUsed'
                },
                {
                    extend: 'excel',
                    filename: 'numSSsnpsUsed',
                    messageTop: 'Number of SNPs Used in Simple Sum Calculation'
                }
            ]
        });
    });
    
}




function list_secondary_SSPvalues(titles, SSPvalues, SSPvaluesN) {
    var tbody = d3.select("#secondary-table").select("tbody");

    // Column headers
    d3.select('#dataset-desc').text('Dataset description');
    secondaryTable.select('thead').select('tr')
        .append('th')
        .attr('class', 'th-sm')
        .text('Simple Sum -log10P');
    secondaryTable.select('thead').select('tr')
        .append('th')
        .attr('class', 'th-sm')
        .text('Number of SNPs Used in SS Calculation')
    
    
    // Table body:
    for(i=0; i<titles.length; i++) { // for each title description
        var row = tbody.append('tr');
        row.append('td').text(titles[i]);
        row.append('td').text(SSPvalues[i]);
        row.append('td').text(SSPvaluesN[i]);
    }

    // Add DataTables functionality:
    secTable = $(document).ready(function () {
        var thesectable = $('#secondary-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'Secondary_datasets_SS_pvalues'
                },
                {
                    extend: 'excel',
                    filename: 'Secondary_datasets_SS_pvalues',
                    messageTop: 'Simple Sum P-values of Secondary Datasets'
                }
                ]
        });
    });
}


