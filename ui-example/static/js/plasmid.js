app.controller('preview', function($scope){
	$scope.name = "";
	$scope.layers = [];
	$scope.markers = [];
	$scope.range = [];
	$scope.width = 25;
	$scope.radius = 100;
	$scope.length = 0;
	$scope.url = '/static/html/plasmid.html'
	$scope.click = function(event, marker){
		part = w2ui['SideBar'].parts.find( ({text}) => ( text == marker.labels[0].text ) );
		w2ui['SideBar'].click(part.dbid);
		// renderPart(part.dbid);
		// // Highlighting selected marker

		// // Dirty hack to remove current legend in sequence-viewer
		// $('.coverageLegend').empty();
		
		// // Getting a number according to dbid last digit - for colouring purposes
		// var number = marker.labels[0].text.split(".")[2]

		// var start = marker.start
		// var end = marker.end
		
		// //var length = (marker.end - marker.start)
		// //var start = (marker.start + length/2)
		// ///var end = start + length
		
		// //console.log("Start and End")
		// //console.log(start, end, length)
		// var feature = [{ start: start , end: end, label: marker.labels[0].text, id: number }];
		// highlight_selected(feature);
	}
	$scope.setPart = function(part){
		$scope.layers = [];
		$scope.range = [0];
		$scope.markers = [];

		$scope.layers.push( [part] );
		$scope.length = part.fullLength;
		$scope.name = part.name;

		partLen = part.length;

		nChildren = part.children.length;
		while (nChildren != 0){
			nChildren = 0;
			layer = [];
			start = 0;
			lastLayer = $scope.layers[ $scope.layers.length-1 ];
			for (i=0; i < lastLayer.length; i++){
				p = lastLayer[i];
				for (j=0; j < p.children.length; j++){
					pp = p.children[j];
					layer.push( pp );
					nChildren ++;
				}
			}
			if (nChildren > 0){
				$scope.layers.push(layer);
				$scope.range.push( $scope.layers.length-1 );
			}
		}

		console.log("Layers: ", $scope.layers);

		for(i=0; i < $scope.layers.length; i++){
			level = [];
			layer = $scope.layers[i];
			start = 0;
			for(j=0; j < layer.length; j++){
				part = layer[j];
				level.push( { start: start - 0.5*partLen, end: start + part.length - 0.5*partLen, text: part.text } );
				start += part.length;
			}
			$scope.markers.push(level);
			console.log(level);
		}
	};
})