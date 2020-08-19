$(document).ready(function() {
    $("#slider-minimalIdentity").slider({
        range: "min",
        min: 0,
        max: 100,
        value: 1,
        slide: function(event, ui) {
            $( "#amountIdentityCutoff" ).val(ui.value );
        },
    });

    $( "#amountIdentityCutoff" ).val($( "#slider-minimalIdentity" ).slider( "value" ) );
});
