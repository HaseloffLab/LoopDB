from loopDB import *
import config

loopDB = LoopDB( config.DATABASE_URL, clean = True  )
loopDB.initFromFile( config.SCHEMA_PATH )
