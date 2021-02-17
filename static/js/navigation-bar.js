var version = "1.4.6 alpha"
var myjumbotron = d3.select("#myjumbotron").append('div').classed('col-md-12 jumbotron text-center', true);
var navbar = d3.select("#navigation-bar").append('div').classed('col-md-12 borderit clearfix', true);


myjumbotron.append('h1').append('a').attr('href','/').text('LocusFocus');
myjumbotron.append('p').append('h4').text("Colocalization Testing Across Datasets");
myjumbotron.append('p').append('h5').text(`Version release: ${version}`);


function addLink(thelink, thetext) {
    navbar
        .append('a')
        .attr("href", `${thelink}`)
        .text(`${thetext}`)
        .attr('role','button')
        .classed('btn btn-primary btn-lg', true);
}
addLink("/previous_session", "Session ID");
addLink("https://locusfocus.readthedocs.io", "Documentation");
addLink("/session_id/0d074aff-76f4-42ab-87ef-5e2fa34b04e8", "Example Output");
addLink("mailto:naim.panjwani@sickkids.ca?subject=LocusFocus", "Contact Us");
addLink("https://mailchi.mp/752ab1c4d516/locusfocus", "Subscribe");

var citationButton = navbar.append('button')
    .classed('btn btn-primary btn-lg', true)
    .attr('data-toggle','modal')
    .attr('data-target','#citeModal')
    .text('Citation');
    
var modal = navbar
    .append('div')
    .attr('id', 'citeModal')
    .classed('modal fade', true)
    .attr('role', 'dialog');
var modalContent = modal
    .append('div')
    .classed('modal-dialog modal-lg modal-dialog-centered', true)
        .append('div')
        .classed('modal-content', true);
var modalHeader = modalContent
    .append('div')
    .classed('modal-header',true);
modalHeader
    .append('h4')
    .classed('modal-title',true)
    .text("LocusFocus: Web-based colocalization for the annotation and functional follow-up of GWAS")
modalHeader
    .append('button')
    .attr('type',"button")
    .attr('class',"close")
    .attr('data-dismiss',"modal")
    .html("&times;");
var modalBody = modalContent
    .append('div')
    .classed('modal-body',true);
modalBody
    .append("p")
    .text("Authors:");
modalBody
    .append('p')
    .text("Naim Panjwani, Fan Wang, Scott Mastromatteo, Allen Bao, Cheng Wang, Gengming He, Jiafen Gong, Johanna M. Rommens, Lei Sun, Lisa J. Strug");
var modalFooter = modalContent
    .append('div')
    .classed('modal-footer',true);
modalFooter
    .append('p')
    .append('a')
    .attr('href','https://doi.org/10.1371/journal.pcbi.1008336')
    .text("PLoS Comput Biol 2020 16(10): e1008336.");
modalFooter
    .append('button')
    .attr('type','button')
    .classed("btn btn-default", true)
    .attr('data-dismiss','modal')
    .text("Close");
