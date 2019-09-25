var theTable = d3.select("#variants-table");

function buildTable(genes, tissues, SSPvalues) {
    // console.log(data);
    var tbody = d3.select("#variants-table").select("tbody");

    // Clear table:
    // if ( $.fn.dataTable.isDataTable( '#variants-table' ) ) {
    //     $('#variants-table').DataTable().destroy();
    // }
    tbody.text("")

    // Build column headers:
    for(i=0; i<tissues.length; i++) {
        theTable.select('thead').select('tr')
            .append('th')
            .attr('class', 'th-sm')
            .text(tissues[i]);
    }

    // Add table body:
    for(i=0; i<genes.length; i++) { // for each tissue
        var row = tbody.append('tr');
        row.append('td').text(genes[i]);
        for(j=0; j<tissues.length; j++) { // for each column
            row.append('td').text(SSPvalues[j][i]);
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
        // $('.dataTables_length').addClass('bs-select');
        // new $.fn.dataTable.Buttons( thedatatable, {
        //     buttons: [
        //         'copy', 'csv', 'excel', 'pdf'
        //     ]
        // });
        // thedatatable.buttons().container()
        //     .appendTo( $('.col-sm-6:eq(0)', thedatatable.table().container() ) );
    });
    
}

