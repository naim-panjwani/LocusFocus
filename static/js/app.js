const startingChr = 1;
const startingPos = 205500000;
const endingPos = 206000000;
const genomicWindowLimit = 2e6;
var submitButton = d3.select("#submit-btn");
var errorDiv = d3.select("#error-messages");
var theTable = d3.select("#variants-table");
var gtexTissuesMsgDiv = d3.select("#tissue-select");
var coordinate = "hg19" // default
var gtex_version = "v7"; // default
var gtexurl = `/gtex/${gtex_version}/tissues_list`;
var markerColDiv = d3.select('#snp');
var chromColDiv = d3.select('#chrom');
var posColDiv = d3.select("#pos");
var refColDiv = d3.select("#ref");
var altColDiv = d3.select("#alt");

var locText = d3.select("#locusText").text();
d3.select("#locusText").text(`${locText} (max: ${genomicWindowLimit/1e6} Mbp):`);
d3.select("#locus").attr('placeholder', `${startingChr}:${startingPos}-${endingPos}`);

// FUNCTIONS

// Coordinate system selection change
function coordinateChange(newCoordinate) {
    $('#LD-populations').multiselect('destroy');
    $('#GTEx-tissues').multiselect('destroy');
    $('#region-genes').multiselect('destroy');
    d3.select("#locus").property("value", "1:205500000-206000000");
    if(newCoordinate === "hg38") {
        gtex_version = "v8";
        gtexurl = `/gtex/${gtex_version}/tissues_list`;
        coordinate = "hg38";
        gtexTissuesMsgDiv.text("Select GTEx (V8) Tissues");
        d3.select("#genes-select").text("Select Genes (enter coordinates above to populate)");
        d3.select("#region-genes").text("");
        init();
    }
    else if(newCoordinate.toLowerCase() == "hg19") {
        gtex_version = "v7";
        gtexurl = `/gtex/${gtex_version}/tissues_list`;
        coordinate = "hg19";
        gtexTissuesMsgDiv.text("Select GTEx (V7) Tissues");
        d3.select("#genes-select").text("Select Genes (enter coordinates above to populate)");
        d3.select("#region-genes").text("");
        init();
    }
}

// helper functions:
function askChromInput(chromColDiv) {
    chromColDiv.append("label")
        .attr('for', 'chrom')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Header text corresponding to the chromosome column in your txt/tsv file (primary dataset)")
        .text("Chromosome Column Name:");
    chromColDiv.append("input")
        .attr('class', 'form-control')
        .attr('name', 'chrom-col')
        .attr('type', 'text')
        .attr('onfocus', 'this.value=""')
        .attr('value', '#CHROM')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the chromosome column in your txt/tsv file (primary dataset)");
}
function askPosInput(posColDiv) {
    posColDiv.append('label')
        .attr('for','position')
        .attr('data-toggle','tooltip')
        .attr('title',"Header text corresponding to the basepair coordinate position column in your txt/tsv file (primary dataset)")
        .text("Position Column Name:");
    posColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'pos-col')
        .attr('type','text')
        .attr('onfocus', 'this.value=""')
        .attr('value','BP')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the basepair coordinate position column in your txt/tsv file (primary dataset)");
}
function askRefInput(refColDiv) {
    refColDiv.append('label')
        .attr('for','ref')
        .attr('data-toggle','tooltip')
        .attr('title',"Header text corresponding to the reference allele column in your txt/tsv file (primary dataset)")
        .text("Reference Allele Column Name:");
    refColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'ref-col')
        .attr('type','text')
        .attr('onfocus', 'this.value=""')
        .attr('value','REF')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the reference allele column in your txt/tsv file (primary dataset)");
}
function askAltInput(altColDiv) {
    altColDiv.append('label')
        .attr('for','alt')
        .attr('data-toggle','tooltip')
        .attr('title',"Header text corresponding to the alternate allele column in your txt/tsv file (primary dataset)")
        .text("Alternate Allele Column Name:");
    altColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'alt-col')
        .attr('type','text')
        .attr('onfocus', 'this.value=""')
        .attr('value','ALT')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the alternate allele column in your txt/tsv file (primary dataset)");
}
function askSNPInput(markerColDiv) {
    markerColDiv.append('label')
        .attr('for','markerColname')
        .attr('data-html','true')
        .attr('data-toggle', 'tooltip')
        .attr('title', "<p>Header text corresponding to the variant name column in your txt/tsv file (primary dataset).</p>")
        .text("Marker Column Name:")
    markerColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'snp-col')
        .attr('type', 'text')
        .attr('onfocus', "this.value='")
        .attr('value', "ID")
        .attr('data-toggle', 'tooltip')
        .attr('data-html',"true")
        .attr('title', "<p>Enter the header text corresponding to the variant name column in your txt/tsv file (primary dataset).</p><p>Accepted formats: rs7512462, 1_205899595_T_C_b37</p>")
}


function markerFormatChange(noSNPID) {
    // clear no-marker input fields:
    markerColDiv.html("");
    chromColDiv.html("");
    posColDiv.html("");
    refColDiv.html("");
    altColDiv.html("");

    if(noSNPID.checked) {
        askChromInput(chromColDiv);
        askPosInput(posColDiv);
        askRefInput(refColDiv);
        askAltInput(altColDiv);
    }
    else {
        askSNPInput(markerColDiv);
    }

    // re-initialize popover and tooltip
    $(function () {
        $('[data-toggle="popover"]').popover()
    });
    $(document).ready(function(){
        $('[data-toggle="tooltip"]').tooltip();
    });
}






function init() {
    String.prototype.replaceAll = function(search, replacement) {
        var target = this;
        return target.split(search).join(replacement);
    };
    
    noSNPID = d3.select("#markerCheckbox");
    if(noSNPID.property("checked")) {
        askChromInput(chromColDiv);
        askPosInput(posColDiv);
        askRefInput(refColDiv);
        askAltInput(altColDiv);
    }
    else {
        askSNPInput(markerColDiv);
    }

    
    // Build LD population selections depending on coordinate build chosen:
    var lddiv = d3.select("#LD-populations");
    if(coordinate.toLowerCase() === "hg19") {
        popcodes = ['EUR', 'AFR', 'AMR','ASN']
        popdesc = ['1000 Genomes 2012 EUR', '1000 Genomes 2012 AFR', '1000 Genomes 2012 AMR', '1000 Genomes 2012 ASN']
        lddiv.text('');
        for(var i = 0; i < popcodes.length; i++) {
            lddiv
                .append("option")
                .text(popdesc[i])
                .property('value', popcodes[i])
        }
        lddiv.property("selectedIndex", 0); // EUR selected by default
    }
    else if(coordinate.toLowerCase() === "hg38") {
        popcodes = ['EUR','AFR','AMR','EAS','SAS', 'NFE']
        popdesc = ['1000 Genomes 2018 EUR', '1000 Genomes 2018 AFR', '1000 Genomes 2018 AMR', '1000 Genomes 2018 EAS', '1000 Genomes 2018 SAS', '1000 Genomes 2018 NFE']
        lddiv.text('');
        for(var i = 0; i < popcodes.length; i++) {
            lddiv
                .append("option")
                .text(popdesc[i])
                .property('value', popcodes[i])
        }
        lddiv.property("selectedIndex", 0); // EUR selected by default
    }


    // Auto-complete for genes field:
    // d3.json(`/genenames/${coordinate}`).then(response => {
    //     $( "#gencodeID" ).autocomplete({
    //         source: response
    //     });
    // });
      
    d3.json(gtexurl).then(response => {
        var gtex_tissues = response.map(k => k);

        // Build GTEx tissues multiselect dropdown
        var gtexdiv = d3.select("#GTEx-tissues");
        gtexdiv.text('');
        for(var i = 0; i < gtex_tissues.length; i++) {
            gtexdiv
                .append("option")
                .attr('value', gtex_tissues[i])
                .text(gtex_tissues[i].replaceAll("_"," "));
        }

        // Multi-Select Initialization
        $(document).ready(function() {
        //     $('#LD-populations').multiselect({
        //         enableClickableOptGroups: true,
        //         maxHeight: 400,
        //         buttonWidth: '400px',
        //         checkboxName: function(option) {
        //             return 'multiselect[]';
        //         }
        //     });
            $('#coordinate').multiselect({
                buttonWidth: '200px',
                checkboxName: function(option) {
                    return 'multiselect[]';
                }
            });
            $('#marker-format').multiselect({
                buttonWidth: '400px',
                checkboxName: function(option) {
                    return 'multiselect[]';
                }
            });
            $('#LD-populations').multiselect({
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

    // initialize popover and tooltip
    $(function () {
        $('[data-toggle="popover"]').popover()
    });
    $(document).ready(function(){
        $('[data-toggle="tooltip"]').tooltip();
    });
    
}

init();