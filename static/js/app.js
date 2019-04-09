const startingChr = 1;
const startingPos = 205500000;
const endingPos = 206000000;
const genomicWindowLimit = 2e6;
var submitButton = d3.select("#submit-btn");
var errorDiv = d3.select("#error-messages");
var theTable = d3.select("#variants-table");
var gtex_version = "v7" // default
// gtexurl = `https://gtexportal.org/rest/v1/dataset/tissueSummary?datasetId=gtex_${gtex_version}&format=json`
gtexurl = "http://grch37.rest.ensembl.org/eqtl/tissue/homo_sapiens?content-type=application/json"

var locText = d3.select("#locusText").text();
d3.select("#locusText").text(`${locText} (max: ${genomicWindowLimit/1e6} Mbp):`);
d3.select("#locus").attr('placeholder', `${startingChr}:${startingPos}-${endingPos}`);

// FUNCTIONS
String.prototype.replaceAll = function(search, replacement) {
    var target = this;
    return target.split(search).join(replacement);
};


d3.json('/populations').then(response => {
    var pops = response;
    // console.log(pops);
    d3.json(gtexurl).then(response => {
        // var gtex_tissues = response.tissueSummary; // this code specific to gtexportal.org
        var gtex_tissues = [];
        Object.keys(response).forEach(k => gtex_tissues.push(k));

        // Build populations multiselect dropdown
        for(var i = 0; i < pops['Population Code'].length; i++) {
            var superpop = pops['Super Population Code'][i];
            var superpopSelect = d3.select("#" + superpop);
            if(superpop === "EUR") {
                superpopSelect
                    .append("option")
                    .attr('value', pops['Population Code'][i])
                    .attr('selected','selected')
                    .text("(" + pops['Population Code'][i] + ") " + pops['Population Description'][i]);                
            }
            else {
                superpopSelect
                    .append("option")
                    .attr('value', pops['Population Code'][i])
                    .text("(" + pops['Population Code'][i] + ") " + pops['Population Description'][i]);
            }
        }

        // Build GTEx tissues multiselect dropdown
        var gtexdiv = d3.select("#GTEx-tissues");
        for(var i = 0; i < gtex_tissues.length; i++) {
            gtexdiv
                .append("option")
                // .attr('value', gtex_tissues[i].tissueSiteDetailId) // this code specific to gtexportal.org
                // .text(gtex_tissues[i].tissueSiteDetailId.replaceAll("_"," ")); // this code specific to gtexportal.org
                .attr('value', gtex_tissues[i])
                .text(gtex_tissues[i].replaceAll("_"," "));
        }

        // Multi-Select Initialization
        $(document).ready(function() {
            $('#LD-populations').multiselect({
                enableClickableOptGroups: true,
                maxHeight: 400,
                buttonWidth: '400px',
                checkboxName: function(option) {
                    return 'multiselect[]';
                }
            });
            
            $('#GTEx-tissues').multiselect({
                enableFiltering: true,
                includeSelectAllOption: true,
                maxHeight: 400,
                buttonWidth: '400px',
                checkboxName: function(option) {
                    return 'multiselect[]';
                }
            });

            $('#ld-type').multiselect({
                buttonWidth: '100px',
                checkboxName: function(option) {
                    return 'multiselect[]';
                }
            });
        });
    });
});
