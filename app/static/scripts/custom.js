$(document).ready(function() {
    for (const [key, value] of Object.entries(sliderPercentages)) {
        $("#slider-range" + key).slider({
            range: true,
            min: value["min"],
            max: value["max"],
            values: [value["valueMin"], value["valueMax"]],
            slide: function(event, ui) {
                $("#amount" + key).val(ui.values[0] + " % - " + ui.values[1] + " %");
            },
            change: function(event, ui) {
                make_nucleotide_image();
            }
        });

        $("#amount" + key).val($("#slider-range" + key).slider("values", 0) + " % - " + $("#slider-range" + key).slider("values", 1) + " %");
    }

    $("#slider-forwardReadRange").slider({
        range: true,
        min: 1,
        max: sliderMaxForwardReadRange,
        values: [minRangefw, maxRangefw],
        slide: function(event, ui) {
            $("#amountForwardReadRange").val(ui.values[0] + " - " + ui.values[1]);
        },
        change: function(event, ui) {
            $("#sequenceImageWarning").css("display", "block");
            $("#identityReadsWarning").css("display", "block");
        }
    });

    $("#amountForwardReadRange").val($("#slider-forwardReadRange").slider("values", 0) + " - " + $("#slider-forwardReadRange").slider("values", 1));

    $("#slider-forwardQuality").slider({
        range: "min",
        min: 1,
        max: sliderQualityForwardMaxQuality,
        value: sliderQualityForwardCutoff,
        slide: function(event, ui) {
            $( "#amountQualityForward" ).val(ui.value );
        },
        change: function(event, ui) {
            update_stats();
            $("#sequenceImageWarning").css("display", "block");
            $("#identityReadsWarning").css("display", "block");
        }
    });

    $( "#amountQualityForward" ).val($( "#slider-forwardQuality" ).slider( "value" ) );

    $("#slider-reverseReadRange").slider({
        range: true,
        min: 1,
        max: sliderMaxReverseReadRange,
        values: [minRangerv, maxRangerv],
        slide: function(event, ui) {
            $("#amountReverseReadRange").val(ui.values[0] + " - " + ui.values[1]);
        },
        change: function(event, ui) {
            $("#sequenceImageWarning").css("display", "block");
            $("#identityReadsWarning").css("display", "block");
        }
    });

    $("#amountReverseReadRange").val($("#slider-reverseReadRange").slider("values", 0) + " - " + $("#slider-reverseReadRange").slider("values", 1));

    $("#slider-reverseQuality").slider({
        range: "min",
        min: 1,
        max: sliderQualityReverseMaxQuality,
        value: sliderQualityReverseCutoff,
        slide: function(event, ui) {
            $( "#amountQualityReverse" ).val(ui.value );
        },
        change: function(event, ui) {
            update_stats();
            $("#sequenceImageWarning").css("display", "block");
            $("#identityReadsWarning").css("display", "block");
        }
    });

    $( "#amountQualityReverse" ).val($( "#slider-reverseQuality" ).slider( "value" ) );


    $("#slider-sequenceLengthRange").slider({
        range: true,
        min: sliderMinSequenceLengthRange,
        max: sliderMaxSequenceLengthRange,
        values: [sliderMinSequenceLengthRange, sliderMaxSequenceLengthRange],
        disabled: sliderMinSequenceLengthRange == sliderMaxSequenceLengthRange,
        slide: function(event, ui) {
            $("#amountSequenceLengthRange").val(ui.values[0] + " - " + ui.values[1]);
        },
        change: function(event, ui) {
            make_seq_image();
        }
    });

    $("#amountSequenceLengthRange").val($("#slider-sequenceLengthRange").slider("values", 0) + " - " + $("#slider-sequenceLengthRange").slider("values", 1));
});
