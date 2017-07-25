FormConstructor = function(backbones){
	var self = this;
	this.backbones = backbones;
	
	this.forms = {
		"add:addNew" 		: addNewForm(self.backbones),
		"add:addMultiple"	: uploadFASTAForm(self.backbones),
		"add:assemble" 		: function(){
								return assemblyForm(self.backbones)
							}
	}
}