var tableselect = d3.select('#SSguidance-table');

function buildSSguidanceTable(genes, tissues, SSP, SSP2) {
    var thead = d3.select('#SSguidance-table').select("thead");
    var tbody = d3.select("#SSguidance-table").select("tbody");

    // Clear table:
    if ( $.fn.dataTable.isDataTable( '#SSguidance-table' ) ) {
        var mytable = $('#SSguidance-table').DataTable();
        mytable.destroy();
    }

    thead.text("");
    tbody.text("");

    // Clear table:
    if ( $.fn.dataTable.isDataTable( '#SSguidance-table' ) ) {
        var mytable = $('#SSguidance-table').DataTable();
        mytable.destroy();
    }

    // Column headers
    header = thead.append('tr');
    header
        .append('th')
        .attr('class', 'th-sm')
        .text('Field');
    header
        .append('th')
        .attr('class','th-sm')
        .text('Value');

    var numGTEx = SSP.length;
    var numSecondary = SSP2.length;
    var numNoeQTL = 0;
    var numFirstStage = 0;
    var numTested = 0;
    var numFailed = 0;
    for(i=0; i<tissues.length; i++) { // for each tissue
        for(j=0; j<genes.length; j++) { // for each gene
            if(SSP[i][j] == -1) {
                numNoeQTL += 1
            } else if(SSP[i][j] == -2) {
                numFirstStage += 1
            } else if(SSP[i][j] == -3) {
                numFailed += 1
            } else if(SSP[i][j] > 0) {
                numTested += 1
            } else {
                numFailed += 1
            }
        }
    }
    for(i=0; i<SSP2.length; i++) {
        if(SSP2[i] == -1) {
            numNoeQTL += 1
        } else if(SSP2[i] == -2) {
            numFirstStage += 1
        } else if(SSP2[i] == -3) {
            numFailed += 1
        } else if(SSP2[i] > 0) {
            numTested += 1
        } else {
            numFailed += 1
        }
    }

    var suggested_SSP = -Math.log10(0.05 / numTested);

    // Table body:
    var row = tbody.append('tr');
        row.append('td').text('Total number of secondary datasets (including GTEx)');
        row.append('td').text(numGTEx+numSecondary);
    var row = tbody.append('tr');
        row.append('td').text('Total number of GTEx datasets');
        row.append('td').text(numGTEx);
    var row = tbody.append('tr');
        row.append('td').text('Total number of user-uploaded secondary datasets');
        row.append('td').text(numSecondary);
    var row = tbody.append('tr');
        row.append('td').text('Number of datasets with no eQTL data (-1)');
        row.append('td').text(numNoeQTL);
    var row = tbody.append('tr');
        row.append('td').text('Number of datasets not passing first stage (-2)');
        row.append('td').text(numFirstStage);
    var row = tbody.append('tr');
        row.append('td').text('Number of datasets with computation error (-3)');
        row.append('td').text(numFailed);
    var row = tbody.append('tr');
        row.append('td').text('Number of datasets tested for colocalization');
        row.append('td').text(numTested);
    var row = tbody.append('tr');
        row.append('td').text('Suggested Simple Sum colocalization threshold at alpha 0.05 (-log10P)');
        row.append('td').text(suggested_SSP);

    
    // Add DataTables functionality:
    SSguidanceTable = $(document).ready(function () {
        var thedatatable3 = $('#SSguidance-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'SSguidanceTable'
                },
                {
                    extend: 'excel',
                    filename: 'SSguidanceTable',
                    messageTop: 'Guidance table SS computation'
                }
            ]
        });
    });
}

