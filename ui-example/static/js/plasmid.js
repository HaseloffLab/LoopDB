app.controller('preview', function($scope){
	$scope.name = "";
	$scope.layers = [];
	$scope.markers = [];
	$scope.range = [];
	$scope.width = 20;
	$scope.radius = 100;
	$scope.length = 0;
	$scope.url = '/static/html/plasmid.html'

	$scope.setPart = function(part){
		console.log("RenderPart: ", part);
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
				level.push( { start: start - 0.5*partLen, end: start + part.length - 0.5*partLen, text: part.dbid } );
				start += part.length;
			}
			$scope.markers.push(level);
			console.log(level);
		}
	};
})