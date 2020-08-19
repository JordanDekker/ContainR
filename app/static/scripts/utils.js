var log_data_json = "";
var serverUrl = "http://" + window.location.hostname + ":" + window.location.port;

function download_tsv() {
    disable_buttons("#download_tsv");

    window.location.href = serverUrl + "/api/export_tsv"

    setTimeout(function() {
        enable_buttons("#download_tsv");
    }, 5000);
}

function checkInp(input_list) {
    input_list.forEach(function (s) {
        if (isNaN(x)) {
            alert("Must input numbers");
            return false;
        }
    })
    return true;
}
function clear_errors(){
    $( ".image_error" ).html("");
    $("#BLAST_succes").html("");
    $("#upload_stuff").html("");
    $("#upload_error").html("");
    $("#upload_succes").html("");
}

function disable_buttons(button_id) {
    $(button_id).prop('disabled', true);
}

function enable_buttons(button_id) {
    $(button_id).prop('disabled', false);
}

function disable_all_buttons() {
    $(function(){
        $( ".action_button" ).prop("disabled", true);
    });
}

function enable_all_buttons() {
    $(function(){
        $( ".action_button" ).prop("disabled", false);
    });
}


function make_seq_image(){
    clear_errors();
    disable_all_buttons();
    var min_seq_len = $( "#slider-sequenceLengthRange" ).slider( "values", 0 )
    var max_seq_len = $( "#slider-sequenceLengthRange" ).slider( "values", 1 )

    // if (checkInp([min_seq_len, max_seq_len])){
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/sequenceLength",
        data: {"min_seq_len" : min_seq_len, "max_seq_len" : max_seq_len},
        statusCode: {
            400: function(){
                $("#sequence_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                remove_graphs();
            },
            404: function(){
                $("#sequence_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                $("#sequence_error").html("Internal server error, please contact the developers");

                document.getElementById("sequenceImage").innerHTML = '';
                document.getElementById("sequenceImage").style.height = null;

            }
        },
        success: function( response ) {
            $( "#sequenceImage" ).css("height","370px");
            visualizeSequenceLength(JSON.parse(response));
            update_stats();
        }
    });
}

function make_paired_image(){
    clear_errors();
    disable_all_buttons();
    var filter_paired = $( "#checkPaired").is(":checked");
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/paired",
        data: {"FilterPaired": filter_paired},
        statusCode: {
            400: function(){
                $("#paired_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                remove_graphs();
            },
            404: function(){
                $("#paired_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                $("#paired_error").html("Internal server error, please contact the developers");
                document.getElementById("pairsImage").innerHTML = '';
                document.getElementById("pairsImage").style.height = null;
            }
        },
        success: function(response){
            $( "#pairsImage" ).css("height","370px");
            visualizePairedReads(JSON.parse(response))
        }
    });
}

function make_nucleotide_image(){
    clear_errors();
    disable_all_buttons();
    $("#load_ratios").css("display", "block");
    $("#sequenceImageWarning").css("display", "block");
    $("#identityReadsWarning").css("display", "block");
    var min_A_value = $( "#slider-rangeA" ).slider( "values", 0 )
    var min_T_value = $( "#slider-rangeT" ).slider( "values", 0 )
    var min_G_value = $( "#slider-rangeG" ).slider( "values", 0 )
    var min_C_value = $( "#slider-rangeC" ).slider( "values", 0 )
    var max_A_value = $( "#slider-rangeA" ).slider( "values", 1 )
    var max_T_value = $( "#slider-rangeT" ).slider( "values", 1 )
    var max_G_value = $( "#slider-rangeG" ).slider( "values", 1 )
    var max_C_value = $( "#slider-rangeC" ).slider( "values", 1 )
    var bin_size = $( "#nucleotide_bin_size" ).val()

    $.ajax({
        type: "POST",
        url: serverUrl + "/api/nucleotide",
        data: { "minAValue" : min_A_value, "minTValue": min_T_value, "minGValue" : min_G_value, "minCValue": min_C_value,
              "maxAValue" : max_A_value, "maxTValue": max_T_value, "maxGValue" : max_G_value, "maxCValue": max_C_value,
              },
        statusCode: {
            400: function(){
                $("#nucleotide_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                $("#load_ratios").css("display", "none");
                remove_graphs();
            },
            404: function(){
                $("#nucleotide_error").html("Not found, please report this error to the developers");
                $("#load_ratios").css("display", "none");
            },
            411: function(){
                update_stats();
                $("#nucleotide_error").html("No reads found within given nucleotide percentage range");
                $("#load_ratios").css("display", "none");
                document.getElementById("fwNucleotideImage").innerHTML = '';
                document.getElementById("rvcNucleotideImage").innerHTML = '';
            },
            500: function(){
                $("#nucleotide_error").html("Internal server error, please contact the developers");
                $("#load_ratios").css("display", "none");
                document.getElementById("fwNucleotideImage").innerHTML = '';
                document.getElementById("rvcNucleotideImage").innerHTML = '';
            }
        },
        success: function(response){
            $( "#fwNucleotideImage" ).css("height","370px");
            $( "#rvcNucleotideImage" ).css("height","370px");
            var combined = JSON.parse(response);
            var fw_json = combined.fw_json;
            var rvc_json = combined.rvc_json;
            $("#load_ratios").css("display", "none");

            visualizeNucleotidePercentages(fw_json, "fwNucleotideImage", "Forward Clustered Nucleotide Percentage");
            visualizeNucleotidePercentages(rvc_json, "rvcNucleotideImage", "Reverse Clustered Nucleotide Percentage");
            update_stats();
        }

    });
}

function make_sequence_image(){
    clear_errors();
    disable_all_buttons();
    $("#load_trimming").css("display", "block");
    $("#sequenceImageWarning").css("display", "none");
    var forward_read = $("#toggle-forwardRead").prop("checked")
    var min_read_range = $("#slider-" + (forward_read ? "forward" : "reverse") + "ReadRange").slider("values", 0)
    var max_read_range = $("#slider-" + (forward_read ? "forward" : "reverse") + "ReadRange").slider("values", 1)
    var quality_cutoff = $("#slider-" + (forward_read ? "forward" : "reverse") + "Quality").slider("option", "value")

    $.ajax({
        type: "POST",
        url: serverUrl + "/api/nucleotidePercentage",
        data: {"minReadRange" : min_read_range, "maxReadRange": max_read_range, "forwardRead": forward_read, "qualityCutoff": quality_cutoff},
        statusCode: {
            400: function() {
                $("#sequence_trim_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                $("#load_trimming").css("display", "none");
                remove_graphs();
            },
            404: function() {
                $("#sequence_trim_error").html("Not found, please report this error to the developers");
                $("#load_trimming").css("display", "none");
            },
            500: function() {
                $("#sequence_trim_error").html("Internal server error, please contact the developers");
                $("#load_trimming").css("display", "none");
                document.getElementById("sequenceImage").innerHTML = "";
            }
        },
        success: function(response){
            $("#sequenceTrimImage").css("height", "370px");
            $("#sequenceQualityImage").css("height", "370px");

            $("#load_trimming").css("display", "none");
            visualizeNucleotideAndQualityPercentages(response);
            update_stats();

        }
    });
}

function update_stats(){
    clear_errors();
    // disable_all_buttons();

    var fw_quality_cutoff = $( "#slider-forwardQuality" ).slider("option", "value")
    var rv_quality_cutoff = $( "#slider-reverseQuality" ).slider("option", "value")

    var min_A_value = $( "#slider-rangeA" ).slider( "values", 0 )
    var min_T_value = $( "#slider-rangeT" ).slider( "values", 0 )
    var min_G_value = $( "#slider-rangeG" ).slider( "values", 0 )
    var min_C_value = $( "#slider-rangeC" ).slider( "values", 0 )
    var max_A_value = $( "#slider-rangeA" ).slider( "values", 1 )
    var max_T_value = $( "#slider-rangeT" ).slider( "values", 1 )
    var max_G_value = $( "#slider-rangeG" ).slider( "values", 1 )
    var max_C_value = $( "#slider-rangeC" ).slider( "values", 1 )

    var min_seq_len = $( "#slider-sequenceLengthRange" ).slider( "values", 0 )
    var max_seq_len = $( "#slider-sequenceLengthRange" ).slider( "values", 1 )

        $.ajax({
        type: "POST",
        url: serverUrl + "/api/update_stats",
        data: { "fwQCut" : fw_quality_cutoff, "rvQCut" : rv_quality_cutoff, "minAValue" : min_A_value,
                "minTValue": min_T_value, "minGValue" : min_G_value, "minCValue": min_C_value,
                "maxAValue" : max_A_value, "maxTValue": max_T_value,  "maxGValue" : max_G_value,
                "maxCValue": max_C_value, "minSeqLen" : min_seq_len, "maxSeqLen" : max_seq_len},
        statusCode: {
            400: function(){
                $("#nucleotide_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                remove_graphs();
            },
            404: function(){
                $("#nucleotide_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                $("#nucleotide_error").html("Internal server error, please contact the developers");
            }
        },
        success: function(response){
            enable_all_buttons();
            var obj = JSON.parse(response);
            document.getElementById('unflaggedReadsLabel').innerHTML = obj.unflaggedReads;
            document.getElementById('fwReadsLabel').innerHTML = obj.fwReads;
            document.getElementById('rvReadsLabel').innerHTML = obj.rvReads;
            document.getElementById('avgQualityLabelFw').innerHTML = obj.avgQualityFw
            document.getElementById('avgQualityLabelRv').innerHTML = obj.avgQualityRv
            document.getElementById('avgLengthLabelFw').innerHTML = +(Math.round(obj.avgLengthFw + "e+2") + "e-2");
            document.getElementById('avgLengthLabelRv').innerHTML = +(Math.round(obj.avgLengthRv + "e+2") + "e-2");
            make_paired_image();
        }
    });
}


function make_identity_image() {
    clear_errors();
    disable_all_buttons();
    $("#load_identity").show();
    $("#identityReadsWarning").css("display", "none");
    var pairedReadMatch = $("#paired_read_match").val();
    var pairedReadMismatch = $("#paired_read_mismatch").val();
    var pairedReadOpeningGap = $("#paired_read_opening_gap").val();
    var pairedReadExtendGap = $("#paired_read_extend_gap").val();
    var pairedReadMinimumPercentage = $("#paired_read_minimum_percentage").val();
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/identity",
        data: { "paired_read_match": pairedReadMatch, "paired_read_mismatch": pairedReadMismatch,
            "paired_read_opening_gap": pairedReadOpeningGap, "paired_read_extend_gap": pairedReadExtendGap,
            "paired_read_minimum_percentage": pairedReadMinimumPercentage},
        statusCode: {
            400: function(){
                console.log("400 error")
                $("#identity_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                remove_graphs();
            },
            404: function(){
                console.log("404 error")
                $("#identity_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                console.log("500 error")
                $("#identity_error").html("Internal server error, please contact the developers");
                $("#identityImage").css("height", "0px");
                document.getElementById("identityImage").innerHTML = "";
            }
        },
        success: function(response){
            $("#identityImage").css("height", "370px");
            forwardReverseCompare(JSON.parse(response));
            update_stats();
            $("#load_identity").hide();
        }
    });
}

/*
this gets called when the do_blast button is clicked. it takes the data from 2 checkboxes which are either true or false
do_blast_check if true the blast function is executed.
*/
function do_BLAST(){
    clear_errors();
    disable_buttons("#BLAST_button");
    $("#blast_load").css("display", "block");
    var use_options = $("#useOptionsCB").is(":checked");
    var min_identity = $("#amountIdentityCutoff").val();
    var reward_penalty = $("#reward_penalty option:selected").val().split("-");
    var gap_costs = $("#gap_costs option:selected").val().split("-");

    update_blast_duration(use_options);

    $.ajax({
        type: "POST",
        url: serverUrl + "/api/BLAST",
        data: {"minIdent": min_identity,"reward": reward_penalty[0], "penalty": reward_penalty[1], "open": gap_costs[0], "extend": gap_costs[1], "use_options": use_options},
        statusCode: {
            400: function(){
                $("#BLAST_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                $("#blast_load").css("display", "none");
                enable_buttons("#BLAST_button");
            },
            404: function(){
                $("#BLAST_error").html("Not found, please report this error to the developers");
                $("#blast_load").css("display", "none");
                enable_buttons("#BLAST_button");
            },
            500: function(){
                $("#BLAST_error").html("BLAST+ package not installed or not working correctly");
                $("#blast_load").css("display", "none");
                enable_buttons("#BLAST_button");
            }
        },
        success: function(response){
            $( "BLAST").css("heigth", "370px");
            $("#blast_load").css("display", "none");
            enable_buttons("#download_xml_button");
            $("#resultsPage").css("display", "inline");
            $("#BLAST_succes").html("Your BLAST has been completed");
            return JSON.parse(response);
        }
    })
}

function update_blast_duration(use_options){
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/blastDuration",
        data: {"use_options": use_options},
        statusCode: {
            400: function(){
                $("#BLAST_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                $("#blast_load").css("display", "none");
                enable_buttons("#BLAST_button");
            },
            404: function(){
                $("#BLAST_error").html("Not found, please report this error to the developers");
                $("#blast_load").css("display", "none");
                enable_buttons("#BLAST_button");
            },
            500: function(){
                $("#BLAST_error").html("BLAST+ package not installed or not working correctly");
                $("#blast_load").css("display", "none");
                enable_buttons("#BLAST_button");
            }
        },
        success: function(response){
            let blast_duration = JSON.parse(response).duration;
            // Start the progress bar. Use 1 as the minimum start duration
            progress(blast_duration, Math.max(1, blast_duration), $("#blastProgressBar"));
        }
    });
}

/*
    Show or hide blast options when checkbox is pressed
*/
function toggle_advanced_blast_options(){
    advdiv = document.getElementById("blastAdvanced")
    document.getElementById("useOptionsCB").addEventListener('change', (event) => {

        console.log("CHANGE")
        if (event.target.checked) {
            advdiv.style.visibility = "visible"
            advdiv.style.height = "auto"
        } else {
            advdiv.style.visibility = "collapse"
            advdiv.style.height = "0px"
        }
    })
}

/*
this gets called when the make_sunburst button is clicked. it does what it says, it generates the sunburst graph.
*/
function make_sunburst_image(){
    clear_errors();

    document.getElementById("main").innerHTML = '<div id="chart"><div id="explanation" style="visibility:hidden"><span id="percentage"></span><br>of sequences are part of this section (click on a bacterium to get extra data on it)</div></div>';

    $.ajax({
        type: "POST",
        url: serverUrl + "/api/sunburst",

        data: {remake : "False"},
        statusCode: {
            400: function(){
                $("#sunburst_error").html("No known records are loaded, please make sure you uploaded your files in this session");
            },
            404: function(){
                $("#sunburst_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                $("#sunburst_error").html("Internal server error, please contact the developers");
            }
        },
        success: function(response){
            var jsonderulo = JSON.parse(response);
            document.getElementById("upload_load").style.display = "none";
            $("#Sunburst").css("height", "750px");
            build_sunburst(jsonderulo['sun_json'])
        }
    });
}

function remake_sunburst_image(){
    clear_errors();
    disable_buttons("#remake_sunburst_button");
    document.getElementById("main").innerHTML = '<div id="chart"><div id="explanation" style="visibility:hidden"><span id="percentage"></span><br>of sequences are part of this section (click on a bacterium to get extra data on it)</div></div>';

    $.ajax({
        type: "POST",
        url: serverUrl + "/api/sunburst",

        data: {remake :"True"},
        statusCode: {
            400: function(){
                $("#sunburst_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                enable_buttons("#remake_sunburst_button");
            },
            404: function(){
                $("#sunburst_error").html("Not found, please report this error to the developers");
                enable_buttons("#remake_sunburst_button");
            },
            500: function(){
                $("#sunburst_error").html("Internal server error, please contact the developers");
                enable_buttons("#remake_sunburst_button");
            }
        },
        success: function(response){
            var json_response = JSON.parse(response);

            document.getElementById("upload_load").style.display = "none";
            $("#Sunburst").css("height", "750px");
            build_sunburst(json_response['sun_json']);
            enable_buttons("#remake_sunburst_button");
        }
    });
}

function delete_node(){
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/delete_node",
        data: { "data_json": JSON.stringify(log_data_json) },
        statusCode: {
            400: function(){
                $("#delete_error").html("No known records are loaded, please make sure you uploaded your files in this session");
            },
            404: function(){
                $("#delete_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                $("#delete_error").html("Internal server error, please contact the developers");
            }
        },
        success: function(response){
            make_sunburst_image();
        }
    });
}

/*
 this gets called when the sunburst graph is clicked. it takes a string which contains the taxonomic data of the clicked
 bacteria. at the succes it makes 3 graphs with an overlay.
 */
function make_bacterial_vizualisation(data_json){
    log_data_json = data_json;
    clear_errors();
    $('#main').css('pointer-events', 'none');
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/bacterial",
        data: { "data_json": JSON.stringify(data_json)},
        statusCode: {
            400: function(){
                $("#bacterial_error").html("No known records are loaded, please make sure you uploaded your files in this session");
            },
            404: function(){
                $("#bacterial_error").html("Not found, please report this error to the developers");
            },
            500: function(){
                $("#bacterial_error").html("Internal server error, please contact the developers");
            }
        },
        success: function(response){
            $( "#chart" ).css("height", "370px");
            $("#overlay").css("display", "block");
            $('#main').css('pointer-events', '');
            make_graphs(response)
            }
    });
}

function remove_graphs(){
    document.getElementById("pairsImage").innerHTML = '';
    document.getElementById("fwNucleotideImage").innerHTML = '';
    document.getElementById("rvcNucleotideImage").innerHTML = '';
    document.getElementById("sequenceImage").innerHTML = '';
    document.getElementById("identityImage").innerHTML = '';
    document.getElementById("pairsImage").style.height = null;
    document.getElementById("fwNucleotideImage").style.height = null;
    document.getElementById("rvcNucleotideImage").style.height = null;
    document.getElementById("sequenceImage").style.height = null;
    document.getElementById("identityImage").style.height = null;
}

function load_blast_results(){
    clear_errors();
    $.ajax({
        type: "POST",
        url: serverUrl + "/api/load_records",
        statusCode: {
            400: function(){
                $("#upload_error").html("No known records are loaded, please make sure you uploaded your files in this session");
                $("#upload_load").css("display", "none");
                enable_buttons("#BLAST_button");
            },
            404: function(){
                $("#upload_error").html("Not found, please report this error to the developers");
                $("#upload_load").css("display", "none");
                enable_buttons("#BLAST_button");
            },
            500: function(){
                $("#upload_error").html("BLAST+ package not installed or not working correctly");
                $("#upload_load").css("display", "none");
                enable_buttons("#BLAST_button");
            }
        },
        success: function(response){
            $("#upload_load").css("display", "none");
            $("#upload_succes").html("The BLAST results are processed.");
            return JSON.parse(response);
        }
    })

}

function progress(timeleft, timetotal, $element) {
    var progressBarWidth = timeleft * $element.width() / timetotal;
    var minutes = Math.floor(timeleft / 60);
    var seconds = timeleft % 60;

    $element.find("div").animate({ width: progressBarWidth }, 500).html(minutes + ":" + (seconds < 10 ? "0" + seconds : seconds));

    if (timeleft > 0) {
        setTimeout(function() {
            progress(timeleft - 1, timetotal, $element);
        }, 1000);
    }
}

function redirectPage(button) {
    window.location.href = button.dataset.url;
}
