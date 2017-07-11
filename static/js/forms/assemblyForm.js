assemblyForm = function(backboneList){
	console.log("ASSEMBLY FORM INIT");
	form = {
		init: function(backboneList){
			this.name  =  "assemblyForm";
			this.header = "Assemble new part";
			this.n =  0;
			this.part  =  {children : [], fullLength: 0, length: 0};
			this.tabs  = [
				{id : "tab0", caption : "Select Backbone"}
			];
			this.fields = [];
			this.page = 0;
			this.actions = {
				clear : {caption: "Clear", onClick : function(){
					w2ui['layout'].content( 'main', assemblyForm(backboneList) );
					renderPart(null);
				}},

				next : {caption : "Proceed", onClick : function(){

					socket.emit('getParts', this.record.selector.site3, function(parts){
						console.log(this.part);

						if (form.n == 0 && parts.length == 0){
							this.onValidate = function(event){
								event.errors.push({ field: this.get("selector"), error : "Can't find parts to proceed assembly." });
							}
						}
						else{
							if (parts.length == 0 && this.record.selector.site3 != this.part.backbone.site5){
								this.onValidate = function(event){
									event.errors.push({ field: this.get("selector"), error: "Can't find parts to proceed assembly." });
								}
							}
							else{
								this.onValidate = function(event){};
							}
						}

						if (this.validate(true).length == 0){
							form.tabs[ form.tabs.length-1 ].caption = this.record.selector.text;
							if (form.n == 0){
								form.part.backbone = this.record.selector;
							}
							else{
								form.part.children.push( this.record.selector );
								form.part.length += this.record.selector.length;
							}
							
							form.part.fullLength += this.record.selector.length;
							
							plasmidScope = angular.element( $("#layout_layout_panel_preview") ).scope();
							plasmidScope.setPart( form.part );
							plasmidScope.$apply();


							form.n++;
							for (i=0; i<form.n; i++){
								form.tabs[i].disabled = true;
							}
							form.page = form.n;

							if (parts.length != 0){
								form.fields[0] = {
									field: "selector", type: "list", required: true,
									options:{
										items: parts
									},
									html:{
										page: form.n,
										caption : "Select part " + form.n
									}
								}
								form.tabs.push( {id: "tab" + form.n, caption: "Select part "+form.n} )
							}
							else{
								form.actions.next.caption = "Finalise";
								form.actions.next.onClick = function(){

									this.onValidate = function(event){
										if (w2ui['sideBar'].parts.find( ({name}) => name == this.record["Part name"] )){
											event.errors.push( {field: this.get("Part name"), error: "Part name already exists"} )
										}
									};

									if (this.validate(true).length == 0){
										form.part.name = this.record["Part name"];
										console.log( form.part )
										socket.emit("submitAssembly", form.part, function(response){
											if (response[0] == "OK"){
												console.log("RESPONSE OK");
												renderPartList();
												renderPart(response[1]);
												w2ui['sideBar'].selected = response[1].dbid;
											}
										});
									}
								}
								form.fields = [
									{ field: "Part name", type: "text", required: true,
										html:{
											page: form.n
										}
									}
								];
								form.tabs.push( { id: "tabF", caption: "Set part name" } );
							}

							w2ui['assemblyForm'].destroy();
							w2ui['layout'].content('main', $().w2form(form) );

						}
						else{
							console.log("INVALID");
						}

					}.bind(this));

				}}
			} // actions
			
			if ('assemblyForm' in w2ui){
				w2ui['assemblyForm'].destroy();
			}

			console.log("Assembly form: ", form);
			backbones = [];
			for(i in backboneList){
				backbones = backbones.concat(backboneList[i].backbones);
			}

			form.fields[0] = {
				field: "selector", type: "list", required: true,
				options:{
					items: backbones.filter( ({site3}) => !!site3 )
				},
				html:{
					caption : "Select backbone"
				}
			}
		} // init
	} // form

	form.init(backboneList);

	return $().w2form(form);
} // function