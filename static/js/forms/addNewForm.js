function addNewForm(backbones){
	return $().w2form({
		name : "addNewForm",
		header : "Add new part",
		backbones: [],
		fields : [
			{ name: 'Part name', type: 'text', required: true,
				html:{
					attr: 'style=width:400px'
				}
			},
			{ name: 'GenBank file', type: 'file', required: false,
				options: {
					max : 1,
					onAdd: function(){
						w2ui['addNewForm'].get("Sequence").disabled = true;
						w2ui['addNewForm'].refresh("Sequence");
					},
					onRemove: function(){
						w2ui['addNewForm'].get("Sequence").disabled = false;
						w2ui['addNewForm'].refresh("Sequence");
					}
				},
				html:{
					attr: 'style=width:400px',
				}
			},
			{ name: 'BaseSeq', type : 'list', required: true,
				options:{
					items: backbones,
					renderDrop: function(item){
						return item.text
					}
				},
				html:{
					attr: 'style=width:200px'
				}
			},

			{ name: 'Backbone', type : 'list', required : true,
				options:{
					items : [],
					renderDrop : function(item){
						return item.text
					},
				},
				html:{
					attr: 'style=width:200px'
				}
			},
			{
				name: 'Sequence', type: 'textarea', required: false,
				html:{
					attr: 'style=width:400px;height:200px;resize:none'
				}
			}
		],
		onChange: function(event){
			if (event.target == "BaseSeq"){
				var bckbField = this.get("Backbone")
				bckbField.options.items = event.value_new.backbones;
				this.set("Backbone", bckbField);
				this.record["BaseSeq"] = event.value_new;
				this.refresh();
			}
			
		},
		onValidate: function(event){
			if ( this.record["GenBank file"] == "" && this.record["Sequence"] == "" ){
				event.errors.push( {field: this.get("GenBank file"), error: "Either GB file or Sequence is required"  } );
				event.errors.push( {field: this.get("Sequence"), error: "Either GB file or Sequence is required"  } );
			}
			if (w2ui['sideBar'].parts.find( ({name}) => name == this.record["Part name"] )){
				event.errors.push( {field: this.get("Part name"), error: "Part name already exists"} );
			}
			if (this.record["Sequence"].match(new RegExp("[^ATCGatcg]", "g"))){
				event.errors.push( {field: this.get("Sequence"), error: "Invalid characters detected"  } );
			}
		},
		actions:{
			"Save" : function(event){
				console.log(this.record);
			
				if (this.validate(true).length == 0){
					socket.emit("addL0", this.record, function(response){
						console.log("RESPONSE", response);
						if (response[0] == "OK"){
							w2ui['layout'].content('main', "<i class='fa fa-check-circle-o fa-5x'></i>");
							
							
							renderPartList(function(){
								w2ui['sideBar'].expandParents(response[1].dbid);
								w2ui['sideBar'].click(response[1].dbid);
							});

						}
						else{
							w2alert(response[1]);
						}
					});
				}
			} 
		}
	})
}
