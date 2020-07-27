
function buildColoc2Table(tabledata) {
    var tableselect = d3.select("#coloc2-table");

    // Column headers
    header = tableselect.append('thead').append('tr');
    header
        .append('th')
        .attr('class', 'th-sm')
        .text('ProbeID');
    header
        .append('th')
        .attr('class','th-sm')
        .text('PP.H4.abf');
    
    
    // Table body:
    tbody = tableselect.append('tbody');
    for(i=0; i<tabledata['ProbeID'].length; i++) {
        var row = tbody.append('tr');
        row.append('td').text(tabledata['ProbeID'][i]);
        row.append('td').text(tabledata['PPH4abf'][i]);
    }

    // Add DataTables functionality:
    coloc2Table = $(document).ready(function () {
        $('#coloc2-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'COLOC2_table'
                },
                {
                    extend: 'excel',
                    filename: 'COLOC2_table',
                    messageTop: 'COLOC2 Posterior Probabilites'
                }
                ]
        });
    });
}

