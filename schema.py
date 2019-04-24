from loopDB import LoopDB # Importing the LoopDB module

loopDB = LoopDB( app.config["DATABASE_URL"] , clean = True)
loopDB.initFromFile('schema.json') # Initialising from the schema file