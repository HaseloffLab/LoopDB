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

toolbar = {
	name: "toolbar",
	items :[
		{ id: 'logo', type: 'html', html: '<img src="static/img/Logo.png" style="height:30px">'},
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
		}
	],

	onClick: function(event) {
		console.log(event);
		if (event.target in formConstructor.forms){
			var form = formConstructor.forms[event.target]; 
			form.clear();
			renderPart(null);
			w2ui['layout'].content( 'main', form );
		}
		
	}
};

layout = $().w2layout({
	name: 'layout',
	panels: [
		{ type: 'left', size: 150, resizable: true, content: sidebar},
		{ type: 'main', toolbar: toolbar},
		{ type: 'preview', size: 400}
	]
});