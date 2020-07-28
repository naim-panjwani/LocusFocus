const startingChr = "1";
const startingPos = "205,500,000";
const endingPos = "206,000,000";
const genomicWindowLimit = 2e6;
var submitButton = d3.select("#submit-btn");
var errorDiv = d3.select("#error-messages");
var theTable = d3.select("#variants-table");
var gtexTissuesMsgDiv = d3.select("#tissue-select");
var coordinate = "hg19" // default
var gtex_version = "v7"; // default
var gtexurl = `/gtex/${gtex_version}/tissues_list`;
var markerColDiv = d3.select('#snp');
var variantInputsDiv = d3.select('#variantInputs');
var chromColDiv = d3.select('#chrom');
var posColDiv = d3.select("#pos");
var refColDiv = d3.select("#ref");
var altColDiv = d3.select("#alt");
var statsDiv = d3.select("#statsDiv");
var statsDiv2 = d3.select("#statsDiv2");


var locText = d3.select("#locusText").text();
d3.select("#locusText").text(`${locText} (max: ${genomicWindowLimit/1e6} Mbp):`);
d3.select("#locus").attr('value', `${startingChr}:${startingPos}-${endingPos}`);

// FUNCTIONS

// Coordinate system selection change
function coordinateChange(newCoordinate) {
    $('#LD-populations').multiselect('destroy');
    $('#GTEx-tissues').multiselect('destroy');
    $('#region-genes').multiselect('destroy');
    d3.select("#locus").property("value", `${startingChr}:${startingPos}-${endingPos}`);
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
    chromColDiv.html("");
    chromColDiv.append("label")
        .attr('for', 'chrom')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Header text corresponding to the chromosome column in your txt/tsv file (primary dataset)")
        .text("Chromosome Column Name:");
    chromColDiv.append("input")
        .attr('class', 'form-control')
        .attr('name', 'chrom-col')
        .attr('type', 'text')
        .attr('value', '#CHROM')
        .attr('onfocus',"this.value=''")
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the chromosome column in your txt/tsv file (primary dataset)");
}
function askPosInput(posColDiv) {
    posColDiv.html("");
    posColDiv.append('label')
        .attr('for','position')
        .attr('data-toggle','tooltip')
        .attr('title',"Header text corresponding to the basepair coordinate position column in your txt/tsv file (primary dataset)")
        .text("Position Column Name:");
    posColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'pos-col')
        .attr('type','text')
        .attr('value','POS')
        .attr('onfocus',"this.value=''")
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the basepair coordinate position column in your txt/tsv file (primary dataset)");
}
function askRefInput(refColDiv) {
    refColDiv.html("");
    refColDiv.append('label')
        .attr('for','ref')
        .attr('data-toggle','tooltip')
        .attr('title',"Header text corresponding to the reference allele column in your txt/tsv file (primary dataset)")
        .text("Reference Allele Column Name:");
    refColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'ref-col')
        .attr('type','text')
        .attr('value','REF')
        .attr('onfocus',"this.value=''")
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the reference allele column in your txt/tsv file (primary dataset)");
}
function askAltInput(altColDiv) {
    altColDiv.html("");
    altColDiv.append('label')
        .attr('for','alt')
        .attr('data-toggle','tooltip')
        .attr('title',"Header text corresponding to the alternate allele column in your txt/tsv file (primary dataset)")
        .text("Alternate Allele Column Name:");
    altColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'alt-col')
        .attr('type','text')
        .attr('value','ALT')
        .attr('onfocus',"this.value=''")
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the alternate allele column in your txt/tsv file (primary dataset)");
}
function askSNPInput(markerColDiv) {
    markerColDiv.html("");
    markerColDiv.append('label')
        .attr('for','markerColname')
        .attr('data-html','true')
        .attr('data-toggle', 'tooltip')
        .attr('title', "<p>Header text corresponding to the variant name column in your txt/tsv file (primary dataset).</p><p>Accepted formats: rs7512462, 1_205899595_T_C_b37.</p>")
        .text("Marker Column Name:");
    markerColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'snp-col')
        .attr('type', 'text')
        .attr('value', "ID")
        .attr('onfocus', "this.value=''")
        .attr('data-toggle', 'tooltip')
        .attr('data-html',"true")
        .attr('title', "<p>Enter the header text corresponding to the variant ID column in your txt/tsv file (primary dataset).</p><p>Accepted formats: rs7512462, 1_205899595_T_C_b37</p>");
}
function addVariantInputs() {
    variantInputsDiv.html("");
    var chromColDiv = variantInputsDiv.append('div').attr('class','col-md-3').attr('id','chrom');
    var posColDiv = variantInputsDiv.append('div').attr('class','col-md-3').attr('id','pos');
    var refColDiv = variantInputsDiv.append('div').attr('class','col-md-3').attr('id','ref');
    var altColDiv = variantInputsDiv.append('div').attr('class','col-md-3').attr('id','alt');
    askChromInput(chromColDiv);
    askPosInput(posColDiv);
    askRefInput(refColDiv);
    askAltInput(altColDiv);
}


function askBetaInput(betaColDiv) {
    betaColDiv.html("");
    betaColDiv.append('label')
        .attr('for', 'betaColname')
        .attr('data-html', 'true')
        .attr('data-toggle', 'tooltip')
        .attr('title', "Enter the header text corresponding to the beta column in your txt/tsv file (primary dataset)")
        .text("Beta Column Name:");
    betaColDiv.append('input')
        .attr('class', 'form-control')
        .attr('name', 'beta-col')
        .attr('type','text')
        .attr('value', "BETA")
        .attr('onfocus',"this.value=''")
        .attr('data-toggle', 'tooltip')
        .attr('data-html','true')
        .attr('title', 'Enter the header text corresponding to the beta column in your txt/tsv file (primary dataset)');
}
function askStdErrInput(stderrColDiv) {
    stderrColDiv.html("");
    stderrColDiv.append("label")
        .attr('for', 'stderrColname')
        .attr('data-html','true')
        .attr('data-toggle', 'tooltip')
        .attr('title','Enter the header text corresponding to the standard error column in your txt/tsv file (primary dataset)')
        .text('Standard Error Column Name:');
    stderrColDiv.append('input')
        .attr('class','form-control')
        .attr('name','stderr-col')
        .attr('type','text')
        .attr('value', "SE")
        .attr('onfocus',"this.value=''")
        .attr('data-toggle','tooltip')
        .attr('data-html','true')
        .attr('title','Enter the header text corresponding to the standard error column in your txt/tsv file (primary dataset)');
}
function askNumSamplesInput(numSamplesDiv) {
    numSamplesDiv.html("");
    numSamplesDiv.append("label")
        .attr('for', 'numSamples')
        .attr('data-html','true')
        .attr('data-toggle', 'tooltip')
        .attr('title','Enter the header text corresponding to the number of samples column in your txt/tsv file (primary dataset)')
        .text('Number of Samples Column Name:');
    numSamplesDiv.append('input')
        .attr('class','form-control')
        .attr('name','numsamples-col')
        .attr('id','numsamples-col')
        .attr('type','text')
        .attr('value', "N")
        .attr('onfocus',"this.value=''")
        .attr('data-toggle','tooltip')
        .attr('data-html','true')
        .attr('title','Enter the header text corresponding to the number of samples column in your txt/tsv file (primary dataset)');
}
function askPvalueInput(pvalueColDiv) {
    pvalueColDiv.html("");
    pvalueColDiv.append("label")
        .attr('for', 'p-value')
        .attr('data-html','true')
        .attr('data-toggle', 'tooltip')
        .attr('title','Header text corresponding to the p-value column in your txt/tsv file (primary dataset)')
        .text('P-value Column Name:');
    pvalueColDiv.append('input')
        .attr('class','form-control')
        .attr('name','pval-col')
        .attr('type','text')
        .attr('value', "P")
        .attr('onfocus',"this.value=''")
        .attr('data-toggle','tooltip')
        .attr('data-html','true')
        .attr('title','Enter the header text corresponding to the p-value column in your txt/tsv file (primary dataset)');
}
function askMafInput(mafColDiv) {
    mafColDiv.html("");
    mafColDiv.append("label")
        .attr('for', 'maf')
        .attr('data-html','true')
        .attr('data-toggle', 'tooltip')
        .attr('title','Header text corresponding to the MAF column in your txt/tsv file (primary dataset)')
        .text('MAF Column Name:');
    mafColDiv.append('input')
        .attr('class','form-control')
        .attr('name','maf-col')
        .attr('type','text')
        .attr('value', "MAF")
        .attr('onfocus',"this.value=''")
        .attr('data-toggle','tooltip')
        .attr('data-html','true')
        .attr('title','Enter the header text corresponding to the MAF column in your txt/tsv file (primary dataset)');
}
function askNumCasesInput(studytype) {
    var numCasesDiv = d3.select('#numCases');
    numCasesDiv.html("");
    if(studytype === 'cc') {
        numCasesDiv.append("label")
            .attr('for', 'numcases')
            .attr('data-html','true')
            .attr('data-toggle', 'tooltip')
            .attr('title','Enter the number of cases in the study')
            .text('Number of Cases:');
        numCasesDiv.append('input')
            .attr('class','form-control')
            .attr('name','numcases')
            .attr('id','numcases')
            .attr('type','text')
            .attr('value', 100)
            .attr('onfocus',"this.value=''")
            .attr('data-toggle','tooltip')
            .attr('data-html','true')
            .attr('title','Enter the number of cases in the study')
            .attr('onchange',"checkNumSamplesInput(this.value)");
        numCasesDiv.append('div').attr('class','input-error').attr('id','numSamplesError-message');
    }
}
function askStudytypeInput(studytypeDiv) {
    studytypeDiv.html("");
    studytypeDiv.append("label")
        .attr('for', 'studytype')
        .attr('data-toggle', 'tooltip')
        .attr('title','Select whether the phenotype is quantitative, or whether  it is a case-control design')
        .text('Select Study Type:');
    studytypeDivSelect = studytypeDiv.append('p')
        .append("select")
        .attr('id','studytype')
        .attr('name','studytype')
        .attr('onchange','askNumCasesInput(this.value)');
    studytypeDivSelect
        .append('option')
        .text('Quantitative')
        .property('value','quant');
    studytypeDivSelect
        .append('option')
        .text('Case-control')
        .property('value', 'cc');
}
function askColocInputs() {
    statsDiv.html("");
    statsDiv2.html("");
    var betaColDiv = statsDiv.append('div').attr('class','col-md-3').attr('id','beta');
    var stderrColDiv = statsDiv.append('div').attr('class','col-md-3').attr('id','stderr');
    var numSamplesDiv = statsDiv.append('div').attr('class','col-md-3').attr('id','numSamples');
    var pvalueColDiv = statsDiv.append('div').attr('class','col-md-3').attr('id','p_value');
    var mafColDiv = statsDiv2.append('div').attr('class','col-md-3').attr('id','maf');
    var studytypeDiv = statsDiv2.append('div').attr('class','col-md-3').attr('id','studytype');
    statsDiv2.append('div').attr('class','col-md-3').attr('id','numCases');
    askBetaInput(betaColDiv);
    askStdErrInput(stderrColDiv);
    askNumSamplesInput(numSamplesDiv);
    askPvalueInput(pvalueColDiv);
    askMafInput(mafColDiv);
    askStudytypeInput(studytypeDiv);
}

function inferVariant(snpbox) {
    if(snpbox.checked) {
        variantInputsDiv.html("");
    }
    else {
        addVariantInputs();
    }
    // re-initialize popover and tooltip
    $(function () {
        $('[data-toggle="popover"]').popover()
    });
    $(document).ready(function(){
        $('[data-toggle="tooltip"]').tooltip();
    });
}

function addColoc2Inputs(colocinput) {
    if(colocinput.checked) {
        askColocInputs();
    }
    else {
        statsDiv.html("");
        statsDiv2.html("");
        var pvalueColDiv = statsDiv.append('div').attr('class','col-md-3').attr('id','p_value');
        askPvalueInput(pvalueColDiv);
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
    
    askSNPInput(markerColDiv);
    askChromInput(chromColDiv);
    askPosInput(posColDiv);
    askRefInput(refColDiv);
    askAltInput(altColDiv);

    statsDiv.html("");
    var pvalueColDiv = statsDiv.append('div').attr('class','col-md-3').attr('id','p_value');
    askPvalueInput(pvalueColDiv);
    
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
        $('[data-toggle="tooltip"]').tooltip({
            delay: { "show": 500, "hide": 100 }
        });
    });
    
}

init();