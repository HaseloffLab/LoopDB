var typeColors = {		'1'		: 'rgb(50, 50, 0)',
						'2'		: 'rgb(50, 0, 0)',
						'3'		: 'rgb(0, 50, 75)',
						'4'		: 'rgb(0, 50, 0)',
						'5'		: 'rgb(0, 50, 75)'
};

var typeColorsSelected = {	'0'		: 'rgb(230, 230, 120)',
							'1'		: 'rgb(255, 255, 50)',
							'2'		: 'rgb(250, 0, 50)',
							'3'		: 'rgb(0, 200, 250)',
							'4'		: 'rgb(0, 200, 100)',
							'5'		: 'rgb(200, 80, 80)',
						  	'6'		: 'rgb(50, 150, 200)',
						  	'7'		: 'rgb(100, 240, 50)',
						  	'8'		: 'rgb(255, 180, 100)',
						  	'9'		: 'rgb(150, 150, 200)'
};



function sort_features(features){
	features.sort(function(a, b) {
    	return a.start - b.start;
	});
}

sequenceCoverage = [];
legend = [];

function highlight(features){
	sequenceCoverage = [];
	legend = [];
	var subSeq = "";
	var previous = 0;
	
	var arrayLength = features.length

	for (var i = 0; i < arrayLength; i++) {

		var start = features[i].start;
		var end = features[i].end;
		var label = features[i].label;

		if (start < previous){
			console.log("Skipping label "+label+" due to being inside another feature");
			continue;
		}

		var index = sequenceCoverage.length % 10;

		previous = end;

		color = typeColorsSelected[index];

		sequenceCoverage.push({
			start:		start,
			end:		end,
			bgcolor:	color,
			color:		"black",
			underscore:	false   
		});
		
		legend.push(
			{name: label, color: color, underscore: false}
		);

		subSeq = subSeq + seq.substring(start, end);
		
		$('.coverageLegend').empty();

		sequence.coverage(sequenceCoverage);
		sequence.addLegend(legend);
	}
}
// Do not know to make distinction from first highlighting to the selected highlighting for reset button to work
function highlight_selected(features){
	var sequenceCoverage = [];
	var legend = [];
	var subSeq = "";
	var previous = 0;
	
	var arrayLength = features.length

	for (var i = 0; i < arrayLength; i++) {

		var start = features[i].start;
		var end = features[i].end;
		var label = features[i].label;
		
		if (start < previous){
			console.log("Skipping label "+label+" due to being inside another feature");
			continue;
		}

		var index = sequenceCoverage.length % 10;

		previous = end;

		color = typeColorsSelected[index];

		sequenceCoverage.push({
			start:		start,
			end:		end,
			bgcolor:	color,
			color:		"black",
			underscore:	false
		});
		
		legend.push(
			{name: label, color: color, underscore: false}
		);

		subSeq = subSeq + seq.substring(start, end);
	
		sequence.coverage(sequenceCoverage)
		sequence.addLegend(legend);


	}
}

function reset(){
	// Dirty hack to remove old Legend. Problems when original features are == 0. I think when you assemble, no features are defined nor retained when storing in database.
	$('.coverageLegend').empty();
	sequence.coverage(sequenceCoverage);
	sequence.addLegend(legend);
	console.log(sequenceCoverage);
};	