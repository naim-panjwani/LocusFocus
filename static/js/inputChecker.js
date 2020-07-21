
function loadGenes(build, region) {
    var chrom = parseInt(region.split(':')[0].toLowerCase().replace('chr',''));
    var startbp = parseInt(region.split(':')[1].split('-')[0].replaceAll(',',''));
    var endbp = parseInt(region.split(':')[1].split('-')[1].replaceAll(',',''));
    var genesdiv = d3.select("#region-genes");
    var genesMsgDiv = d3.select("#genes-select");
    d3.json(`/genenames/${build}/${chrom}/${startbp}/${endbp}`).then(response => {
        genesdiv.text('');
        var genenames = response.map(k => k);
        if(genenames.length === 0) {
            genesMsgDiv.text(`No Genes Found in ${region} (${build})`);
        }
        else {
            genesMsgDiv.text(`Select Genes Found in ${region} (${build})`);
        }
        for(var i = 0; i < genenames.length; i++) {
            genesdiv
                .append("option")
                .attr('value', genenames[i])
                .text(genenames[i]);
        }
        $('#region-genes').multiselect('destroy');
        $(document).ready(function() {
            $("#region-genes").multiselect({
                enableFiltering: true,
                includeSelectAllOption: true,
                maxHeight: 400,
                buttonWidth: '400px',
                checkboxName: function(option) {
                    return 'multiselect[]';
                }
            });
        });
    });
}

function checkLocusInput(regiontext) {
    var build = d3.select("#coordinate").property("value");
    var errorSelect = d3.select("#locusErrorDiv");
    errorSelect.text("");
    if(regiontext !== "") {
        d3.json(`/regionCheck/${build}/${regiontext}`).then(response => {
            var message = response['response'];
            if(message !== "OK") {
                errorSelect.text(message);
            }
            else {
                // var genesButton = d3.select("#load-genes-button");
                // genesButton.property('value', `Load Genes in ${regiontext}`);
                loadGenes(build, regiontext);
            }
        })
    }
}

function checkSSInput(regiontext) {
    var build = d3.select("#coordinate").property("value");
    var errorSelect = d3.select("#locusErrorDiv");
    errorSelect.text("");
    if(regiontext !== "") {
        d3.json(`/regionCheck/${build}/${regiontext}`).then(response => {
            var message = response['response'];
            if(message !== "OK") {
                errorSelect.text(message);
            }
        })
    }
}
