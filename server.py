from loopDB import *
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, make_response
from flask_socketio import SocketIO
import json
import base64
from cStringIO import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import Restriction
from domesticate import domesticate, partColors
from Bio.SeqFeature import SeqFeature, FeatureLocation

import config

app = Flask(__name__)
app.config.from_object('config')

socketio = SocketIO(app)

loopDB = LoopDB( app.config["DATABASE_URL"] )

def w2uiFormToDict(request):
	return json.loads(request.form.to_dict()['request'])["record"]

def partToJson(part):
	j = {}
	j["children"] 	= [ partToJson(child) for child in part.children ]
	j["dbid"]		= part.dbid
	j["level"]		= part.level
	j["length"]		= len(part)
	j["fullLength"] = len(part.fullSeq)
	j["name"] 		= part.name
	j["text"] 		= part.name if part.name else part.dbid
	j["site5"]		= part.backbone.adapter.site5
	j["site3"]		= part.backbone.adapter.site3
	j["backbone"]   = backboneToJSON(part.backbone)
	j["color"] 		= partColors[part.backbone.adapter.name]
	return j

def softDelete(dbid):
	session = loopDB.Session()
	
	part = session.query(Part).filter(Part.dbid == dbid).first()
	childPartship = session.query(Partship).filter( Partship.childID == part.id ).all()
	
	if childPartship:
		return ["Error", "Dependent parts exist"]
	else:
		parentPartship = session.query(Partship).filter( Partship.parentID == part.id ).all()
		
		for partship in parentPartship:
			session.delete(partship)

		features = session.query(Feature).filter(Feature.partID == part.id).delete()
		
		print "Soft deleting ", part.dbid

		session.delete(part)
		session.commit()
		return ["OK", part.dbid]

def delete(dbid):
	session = loopDB.Session()
	
	part = session.query(Part).filter(Part.dbid == dbid).first()

	childPartship = session.query(Partship).filter( Partship.childID == part.id ).all()
	parentPartship = session.query(Partship).filter( Partship.parentID == part.id ).all()

	for partship in childPartship:
		parent = session.query(Part).filter(Part.id == partship.parentID).first()
		session.delete(partship)
		delete(parent.dbid)

	for partship in parentPartship:
		session.delete(partship)

	features = session.query(Feature).filter(Feature.partID == part.id).delete()
	session.delete(part)
	session.commit()

	return["OK", part.dbid]

def editName(dbid, newName):
	session = loopDB.Session()
	part = session.query(Part).filter(Part.dbid == dbid).first()
	
	if part:
		part.name = newName
		retMessage = ["OK", partToJson(part)]
	else:
		retMessage = ["Error", "Part not found."]
	
	try:
		session.commit()
	except:
		retMessage = ["Error", "Bad name. Try again."]

	session.close()
	return retMessage

def featureToJson(feature):
	j = {}
	if "label" in feature.qualifiers:
		j["label"] 	= feature.qualifiers["label"]
	else:
		j["label"]  = ""

	j["start"]		= int(feature.location.start)
	j["end"]		= int(feature.location.end)
	j["forward"]	= feature.location.strand
	
	if "ApEinfo_fwdcolor" in feature.qualifiers:
		j["color"]	= feature.qualifiers["ApEinfo_fwdcolor"]
	else:
		j["color"]  = ""
		
	return j

def baseSeqToJSON(baseSeq):
	j = {}
	j["name"]		= baseSeq.name
	j["text"]		= baseSeq.name if baseSeq.name else baseSeq.dbid
	j["dbid"]		= baseSeq.dbid
	j["backbones"]	= list( map( backboneToJSON, baseSeq.backbones ) )
	return j

def backboneToJSON(backbone):
	j = {}
	j["name"]		= backbone.name
	j["text"]		= backbone.name if backbone.name else backbone.dbid
	j["dbid"]		= backbone.dbid
	j["length"]    	= len(backbone.seq)
	if backbone.baseSeq.receiver:
		j["site3"]		= backbone.baseSeq.receiver.site5
		j["site5"]		= backbone.baseSeq.receiver.site3
	
	return j

def digest(part, enzyme):
	seq = part.fullSeq
	if hasattr(Restriction, enzyme):
		enzyme = getattr(Restriction, enzyme)
		
@app.route('/')
def index():
	return render_template('layout.html')

@app.route('/export')
def export():
	dbid = request.args.get('dbid','')

	if dbid:
		part = loopDB.session.query(Part).filter(Part.dbid == dbid ).first()

	response = make_response(part.fullRecord.format("gb"))
	response.headers["Content-Type"] = "application/octet-stream"
	response.headers["Content-Disposition"] = "attachement; filename={0}".format(part.name+'.gb')
	return response

@socketio.on('addL0')
def addL0(part):
	backbone = loopDB.session.query(Backbone).filter(Backbone.dbid == part["Backbone"]["dbid"]).first()

	if part["Sequence"]:
		record = SeqRecord( seq = Seq( part["Sequence"], IUPAC.unambiguous_dna ) )
	else:
		gbFile = StringIO(base64.decodestring( part["GenBank file"][0]["content"] ) )
		record = SeqIO.read(gbFile, format="genbank")

	record = Seq( backbone.adapter.site5, IUPAC.unambiguous_dna)\
				+ record + Seq( backbone.adapter.site3, IUPAC.unambiguous_dna)

	record = domesticate(record)
	#for feature in record.features:
		#print feature.id, feature.qualifiers

	if record:
			record = record[len(backbone.adapter.site5):-len(backbone.adapter.site3)]

			record.features.insert(0, SeqFeature( FeatureLocation(0, len(record)),
				type = "misc_feature", id = "Part Feature", strand = 1,
					qualifiers = {"label" : [part["Part name"]], "ApEinfo_fwdcolor": [ partColors[backbone.adapter.name] ] } ) )

			newPart = loopDB.addPart(backbone = backbone, name = part["Part name"], record = record)
			#print "New Part: ", newPart
			loopDB.commit()
			if newPart:
				return ["OK", partToJson(newPart)]
			else:
				return ["Error", "Failed to insert part to LoopDB."]
	else:
		return ["Error", "Automatic domestication failed. Please domesticate your sequence manualy."]

@socketio.on('getParts')
def getParts(site3):

	session = loopDB.Session()
	if site3 == '':
		parts = session.query(Part).all()
	else:
		parts = session.query(Part).filter(Part.backboneID == Backbone.id).\
				filter(Backbone.adapterID == RES.id).filter(RES.site5 == site3).all()
	jParts = list( map(partToJson, parts) )
	session.close()
	return jParts

@socketio.on('getRecord')
def getRecord(dbid):
	session = loopDB.Session()
	part = session.query(Part).filter(Part.dbid == dbid).first()
	
	if part:
		record = {'seq' : str(part.record.seq), 'features' : list( map(featureToJson, part.record.features) )}
	else:
		record = {}

	session.close()

	return record

@socketio.on('getBackbones')
def getBackbones():
	session = loopDB.Session()
	baseSeq = session.query(BaseSeq).all()
	jBackbones = list( map( baseSeqToJSON, baseSeq ) )
	#print jBackbones
	session.close()
	return jBackbones

@socketio.on('submitAssembly')
def submitPart(part):

	childrenIDs = [ ch["dbid"] for ch in part["children"] ]
	children = loopDB.session.query(Part).filter( Part.dbid.in_( childrenIDs )).all()

	children = [ next( ch for ch in children if ch.dbid == dbid ) for dbid in childrenIDs  ]

	backbone = loopDB.session.query(Backbone).filter(Backbone.dbid == part["backbone"]["dbid"]).first()
	if len(children) == len(part["children"]) and backbone:

		newPart = loopDB.addPart(name = part["name"], children = children, backbone = backbone)
		loopDB.commit()
		return ["OK", partToJson(newPart)]
	else:
		return ["ERROR", "Backbone or one of the children not found"]

@socketio.on('deletePartSoft')
def deletePartSoft(dbid):
	return softDelete(dbid)

@socketio.on('deletePartHard')
def deletePartHard(dbid):
	return delete(dbid)

@socketio.on('editName')
def editPartName(dbid, newName):
	return editName(dbid, newName)

if __name__ == '__main__':
	socketio.run(app, debug=True,host='0.0.0.0', port = 8000)
