import itertools

from partsdb import PartsDB
from tables import *

from sqlalchemy.schema import DropTable
from sqlalchemy.ext.compiler import compiles
from sqlalchemy import or_

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import json

@compiles(DropTable, "postgresql")
def _compile_drop_table(element, compiler, **kwargs):
    return compiler.visit_drop_table(element) + " CASCADE"

class LoopDB(PartsDB):
	def __init__(self, address, clean = False):
		super(LoopDB, self).__init__(address, Base = Base, clean = clean)

	def initFromFile(self, fileName):
		file = open(fileName)
		data  = json.loads( file.read() )

		reDict = {}

		for re in data["RE"]:
			reDict[re["name"]] = self.addRE( **re )

		resDict = {}
		for res in data["RES"]:
			res["re"] = reDict[res["re"]]
			resDict[res["name"]] = self.addRES(**res)

		baseSeqDict = {}
		for baseSeq in data["BaseSeq"]:
			if "gbFile" in baseSeq:
				baseSeq["gbFile"] = open( baseSeq["gbFile"] )
			if "receiver" in baseSeq:
				baseSeq["receiver"] = resDict[ baseSeq["receiver"] ]
			baseSeqDict[ baseSeq["name"] ] = self.addBaseSeq(**baseSeq)

		for backbone in data["Backbone"]:
			backbone["adapter"] = resDict[ backbone["adapter"] ]
			backbone["baseSeq"] = baseSeqDict[ backbone["baseSeq"] ]
			self.addBackbone(**backbone)
		self.commit()

	def _resTest(self, part):
		if part.children:
			reses = [ part.receiverSites[0] ] + \
						list(itertools.chain.from_iterable( \
							[ ch.sites for ch in part.children ] )) + \
								[ part.receiverSites[1] ]

			return not False in [ reses[2*i] == reses[2*i + 1] for i in range( len(reses) / 2 ) ]

		return True

	partTests = [ lambda self, part: part, _resTest]

	def verifyPart(self, part):
		for test in self.partTests:
			if not test(self, part):
				raise Exception("Verification of part {0} failed. Pelase chech your overhangs.".format(part.name))
		return True

	def extractFeatures(self, record):
		dbFeatures = []
		
		for feature in record.features:
			dbFeature = self.addFeature(  start = feature.location.start,
											end = feature.location.end,
											forward = True if feature.location.strand == 1 else False,
											type = feature.type,
											label = feature.qualifiers["label"][0] if "label" in feature.qualifiers else\
														(feature.qualifiers["ApEinfo_label"][0] if "ApEinfo_label" in feature.qualifiers else feature.id),
											color = feature.qualifiers["ApEinfo_fwdcolor"][0] if "ApEinfo_fwdcolor" in feature.qualifiers else "")
			dbFeatures.append(dbFeature)
		return dbFeatures

	def extractBaseFeatures(self, record):
		dbFeatures = []
		
		for feature in record.features:
			dbFeature = self.addBaseFeature(  start = feature.location.start,
											end = feature.location.end,
											forward = True if feature.location.strand == 1 else False,
											type = feature.type,
											label = feature.qualifiers["label"][0] if "label" in feature.qualifiers else\
														(feature.qualifiers["ApEinfo_label"][0] if "ApEinfo_label" in feature.qualifiers else feature.id),
											color = feature.qualifiers["ApEinfo_fwdcolor"][0] if "ApEinfo_fwdcolor" in feature.qualifiers else "")
			dbFeatures.append(dbFeature)
		return dbFeatures

	def _get(self, Table, name):
		entry = self.session.query(Table).filter( or_(Table.name == name,Table.dbid == name) ).first()
		if not entry:
			raise NameError("{0} with name '{1}' not found".format( Table.__tablename__, name ))
		else:
			return entry

	def _tryToUpdate(self, Table, **kwargs):
		row = self.session.query(Table).filter(Table.name == kwargs["name"]).first()
		if row:
			for key, value in kwargs.iteritems():
				setattr(row, key, value)
			self.session.commit()
		
		return row
	
	def getRE(self, name):
		return self._get(RE, name)

	def getRES(self, name):
		return self._get(RE, name)

	def getBaseSeq(self, name):
		return self._get(BaseSeq, name)

	def getBackbone(self, name):
		return self._get(Backbone, name)

	def getPart(self, name):
		return self._get(Part, name)

	def addRE(self, *args, **kwargs):
		re = self._tryToUpdate(RE, **kwargs)
		
		if not re:
			re = super(LoopDB, self).addPart('re', *args, **kwargs)
		
		return re

	def addRES(self, *args, **kwargs):
		
		if "re" in kwargs and isinstance(kwargs["re"], str):
			kwargs["re"] = self._get(RE, kwargs["re"])

		res = self._tryToUpdate(RES, **kwargs)
		
		if not res:
			res = super(LoopDB, self).addPart('res', *args, **kwargs)
		
		return res

	def addBaseSeq(self, *args, **kwargs):
		
		if "gbFile" in kwargs:
			gbFile = kwargs["gbFile"]
			kwargs.pop("gbFile", None)
			record = SeqIO.read(gbFile, format="genbank")
			kwargs["features"] = self.extractBaseFeatures(record)
			kwargs["seq"] = str(record.seq)

		if "receiver" in kwargs and isinstance(kwargs["receiver"], str):
			kwargs["receiver"] = self._get(RES, kwargs["receiver"])

		baseFeature = self.addBaseFeature( start = 0, end = len( kwargs["seq"] ), label = kwargs["name"], type = "misc_feature", forward = True, color = "#FFFF00")
		if "features" in kwargs:
			kwargs["features"].append(baseFeature)
		else:
			kwargs["features"] = [ baseFeature ]

		baseSeq = self._tryToUpdate(BaseSeq, **kwargs)
		
		if not baseSeq:
			baseSeq = super(LoopDB, self).addPart('baseseq', *args, **kwargs)
		
		return baseSeq 
	
	def addBackbone(self, *args, **kwargs):

		if "baseSeq" in kwargs and isinstance(kwargs["baseSeq"], str):
			kwargs["baseSeq"] = self._get(BaseSeq, kwargs["baseSeq"])

		if "adapter" in kwargs and isinstance(kwargs["adapter"], str):
			kwargs["adapter"] = self._get(RES, kwargs["adapter"])

		backbone = self._tryToUpdate(Backbone, **kwargs)
		
		if not backbone:
			backbone = super(LoopDB, self).addPart('backbone', *args, **kwargs)

		return backbone

	def addPart(self, *args, **kwargs):
		
		if "backbone" in kwargs and isinstance(kwargs["backbone"], str):
			kwargs["backbone"] = self._get(Backbone, kwargs["backbone"])

		children = []
		if "children" in kwargs:
			children = kwargs["children"]
			kwargs.pop("children", None)

		if "record" in kwargs:
			record = kwargs["record"]
			kwargs.pop("record", None)
			kwargs["features"] = self.extractFeatures(record)
			kwargs["seq"] = str(record.seq)
			kwargs["seq"] = kwargs["seq"].upper()


		part = self._tryToUpdate(Part, **kwargs)
		if not part:
			newPart = super(LoopDB, self).addPart('part', *args, **kwargs)
			
			for pos, child in enumerate(children):
				if isinstance(child, str):
					child = self._get(Part, child)
				partship = Partship( parent = newPart, child = child, pos = pos )

			if self.verifyPart( newPart ):
				return newPart
			else:
				return None
		else:
			return part
		

	def addFeature(self, *args, **kwargs):
		return super(LoopDB, self).addPart('feature', *args, **kwargs)
	
	def addBaseFeature(self, *args, **kwargs):
		return super(LoopDB, self).addPart('basefeature', *args, **kwargs)



