// Initialising plasmidPreviewController from plasmid.js
angularPlasmid = angular.module("angularPlasmid", ["angularplasmid"]);
angularPlasmid.controller("plasmidPreview", plasmidPreviewController);

socket = io();

formConstructor = null;
downloadFrame = document.createElement('iframe');

downloadPart = function(){
	part = w2ui['sideBar'].parts.find( ({dbid}) => ( dbid == w2ui["sideBar"].selected ) );
	console.log(part);
	if ( part != undefined ){
		downloadFrame.src = "/export?dbid=" + part.dbid;
	}
}

editPart = function(){
	part = w2ui['sideBar'].parts.find( ({dbid}) => ( dbid == w2ui["sideBar"].selected ) );
	if ( part!= undefined ){
		w2prompt({
			label	: "New Name",
			value	: part.text,
			title	: "Change Part Name",
			ok_text	: "Submit",
			cancel_text	: "Cancel"
		})
		.ok( function(event){
			socket.emit("editName", part.dbid, event, function(response){
				part.text = event;
				part.name = event;
				renderPartList(function(){
					w2ui['sideBar'].expandParents(response[1].dbid);
					w2ui['sideBar'].click(response[1].dbid);
				});
			});
		} );
	}
}

deletePart = function(){
	part = w2ui['sideBar'].parts.find( ({dbid}) => ( dbid == w2ui["sideBar"].selected ) );
	console.log(part);
	w2confirm({
		msg: 'Are you sure you want to delete ' + part.text + '?',
		title: 'DNALooper'
	})
	.yes(function(){
		socket.emit('deletePartSoft', part.dbid, function(response){
			if (response[0] == "OK"){
				w2ui['sideBar'].remove(response[1]);
				renderPart(null);
			}
			else{
				w2confirm({
					msg: 'Deleting this part will also delete some parts that are dependent on it. Do you want to proceed ?',
					title: 'Warning!'
				})
				.yes(function(){
					socket.emit('deletePartHard', part.dbid, function(response){
						w2ui['sideBar'].remove(response[1]);
						renderPart(null);
					})
				})
				.no(function(){});
			}
		});
	})
	.no(function(){});
}

$(function () {

	// Adding download support
	downloadFrame.style.display = 'none';
	document.body.appendChild(downloadFrame);

	// Rendering layout defined in layout.js
	$('#layout').w2render('layout');
	w2ui['layout'].load('main', '/static/html/intro.html');
	w2ui['layout'].load('bottom', '/static/html/footer.html');
	// Setting bindings for angularPlasmid
	$( w2ui['layout'].el('preview') ).attr("ng-include", "url");

	angular.element( function() {
		angular.bootstrap( document, ['angularPlasmid'] );
	});

	// Getting Schema backbones
	socket.emit('getBackbones', function(backbones){
		formConstructor = new FormConstructor(backbones);
		console.log(formConstructor);
		renderPartList();
	});
})