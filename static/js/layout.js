sidebar = $().w2sidebar({
	name: "sideBar",
	nodes: [],
	part: null,
	onClick: function(event){
		if( event.object.part){
			part = this.parts.find( ({dbid}) => ( dbid == event.object.id ) );
			renderPart(part);
		}
	}
});

searchField = '<div id="searchField" class="w2ui-field"><label style="position:relative;top:-6px;left:-3px;width:30px"><i class="fa fa-search" aria-hidden="true"></i></label><input type="combo"></div>'

leftBar = $().w2layout({
	name: "leftBar",
	panels: [
		{type: "top", size: 40, content: searchField, overflow: 'hidden'},
		{type: "main", content: sidebar}
	]
});

toolbar = {
	name: "toolbar",
	items :[
		{ id: 'logo', type: 'html', html: '<img src="static/img/LoopDesigner.png" style="height:30px">'},
		{ id: 'spacer', type: 'spacer'},
		{ id: 'add', type: 'menu', caption: 'Add', img: 'w2ui-icon-plus', 
			items:[
				{
					id: "addNew",
					form: "addNewForm",
					text : "Add part",
					icon: "fa fa-circle-o-notch"
				},
				{
					id: "addMultiple",
					text : "Add multiple parts",
					icon: "fa fa-list"
				},
				{
					id: "assemble",
					text : "Assemble part",
					icon: "fa-infinity"
				}
			]
		},
		{ type: 'break' },
		{ id: 'about', type: 'button', caption: 'About', icon: 'fa fa-info', onClick: function(event){
				w2ui['layout'].load('main', '/static/html/about.html');
			}
		},
		{ id: 'help', type: 'button', caption: 'Help', icon: 'fa fa-question', onClick: function(event){
				w2ui['layout'].load('main', '/static/html/help.html');
			}
		}
	],

	onClick: function(event) {
		console.log(event);
		if (event.target in formConstructor.forms){
			
			var form = null;

			//Dynamic forms
			if (typeof formConstructor.forms[event.target] === 'function'){
				form = formConstructor.forms[event.target]();
			}
			//Static forms
			else{
				form = formConstructor.forms[event.target]; 
			}
			
			console.log("FORM: ", form);

			form.clear();

			renderPart(null);
			w2ui['layout'].content( 'main', form );
		}
		
	}
};

layout = $().w2layout({
	name: 'layout',
	panels: [
		{ type: 'left', size: 300, resizable: true, content: leftBar, minSize:232},
		{ type: 'main', toolbar: toolbar},
		{ type: 'preview', resizable: true, size: 300},
		{ type: 'bottom', size: 56}
	]
});