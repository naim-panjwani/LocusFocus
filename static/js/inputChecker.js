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

function checkNumSamplesInput(numsamples) {
    errordiv = d3.select("#numSamplesError-message");
    errordiv.text("");
    if(Number.isInteger(+numsamples) === false) {
        errordiv.text("Must be integer")
    }
}