function buildParamsTable(data, sessionid) {
    var tableselect = d3.select('#params-table');

    // Column headers
    header = tableselect.append('thead').append('tr');
    header
        .append('th')
        .attr('class', 'th-sm')
        .text('Field');
    header
        .append('th')
        .attr('class','th-sm')
        .text('Value');
    
    
    // Table body:
    tbody = tableselect.append('tbody');
    var row = tbody.append('tr');
        row.append('td').text('Session ID');
        row.append('td').text(sessionid);
    var row = tbody.append('tr');
        row.append('td').text('Lead SNP');
        row.append('td').text(data['lead_snp']);
    var row = tbody.append('tr');
        row.append('td').text('Chromosome');
        row.append('td').text(data['chrom']);
    var row = tbody.append('tr');
        row.append('td').text('Start position');
        row.append('td').text(data['startbp']);
    var row = tbody.append('tr');
        row.append('td').text('End position');
        row.append('td').text(data['endbp']);
    var row = tbody.append('tr');
        row.append('td').text('Build');
        row.append('td').text(data['coordinate']);
    var row = tbody.append('tr');
        row.append('td').text('Infer variants');
        row.append('td').text(data['inferVariant']);
    var row = tbody.append('tr');
        row.append('td').text('Number of SNPs');
        row.append('td').text(data['snps'].length);
    var row = tbody.append('tr');
        row.append('td').text('LD Population');
        row.append('td').text(data['ld_populations']);
    var row = tbody.append('tr');
        row.append('td').text('GTEx version');
        row.append('td').text(data['gtex_version']);
    var row = tbody.append('tr');
        row.append('td').text('Number of GTEx tissues selected');
        row.append('td').text(data['gtex_tissues'].length);
    var row = tbody.append('tr');
        row.append('td').text('Number of GTEx genes selected');
        row.append('td').text(data['gtex_genes'].length);
    var row = tbody.append('tr');
        row.append('td').text('SS region');
        row.append('td').text(data['SS_region']);
    var row = tbody.append('tr');
        row.append('td').text('First stage -log10(SS P-value) threshold');
        row.append('td').text(data['set_based_p']);
    var row = tbody.append('tr');
        row.append('td').text('Many SNPs not matching GTEx SNPs');
        row.append('td').text(data['snp_warning']);
    var row = tbody.append('tr');
        row.append('td').text('SNPs matching threshold level');
        row.append('td').text(data['thresh']);
    var row = tbody.append('tr');
        row.append('td').text('Number of SNPs matching with GTEx');
        row.append('td').text(data['numGTExMatches']);
    var row = tbody.append('tr');
        row.append('td').text('Number of user-provided secondary datasets');
        row.append('td').text(data['secondary_dataset_titles'].length);
    var row = tbody.append('tr');
        row.append('td').text('Run COLOC2');
        row.append('td').text(data['runcoloc2']);
        
    

    

    // Add DataTables functionality:
    paramsTable = $(document).ready(function () {
        $('#params-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'parameters_table'
                },
                {
                    extend: 'excel',
                    filename: 'parameters_table',
                    messageTop: 'Input parameters'
                }
                ]
        });
    });
}

