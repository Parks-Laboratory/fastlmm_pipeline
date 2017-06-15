import pyodbc
import csv
import os
from decimal import Decimal
from functools import reduce
import argparse

# Create parser to receive instructions from command line
parser = argparse.ArgumentParser(description = 'Input arguments to populate database')
# Argument to create table and whether to create new table
parser.add_argument('tablename', action = 'store', help = "Table for inserting data")
parser.add_argument('-c', '--create', action = 'store_true', help = "Create table", default= False) 
# Argument to specify path
parser.add_argument('-p', '--path', action = 'store', help = "Directory path with txt file", default=".")
# Argument to specify database
parser.add_argument('-db', '--database', action = 'store', help="Database to be opened", default = "Epistasis")
# Parse all the the arguments together
args = parser.parse_args()
# Tokenize each argument into variables
tablename = args.tablename
create = args.create
path = args.path
database = args.database

# method to create table if needed
def createTable(database, tablename):
	global cursor
	
	query = "create table {!s} ".format(tablename) + "(" + \
			"Trait varchar(max)," \
			" SNP char(50)," \
			" P float," \
			" Beta float," \
			" CONSTRAINT EID PRIMARY KEY (SNP));"
	print(query)
	
	cursor.execute(query)
	cursor.commit()
	print("table %s successfully created in database %s" %(tablename, database))

# method to connect to the database
def createConnection(server, database):
	cn = pyodbc.connect('DRIVER={SQL Server}' + \
						';SERVER=' + server + \
						';DATABASE=' + database + \
						';Trusted_Connection= Yes')
	return cn

# "main" method which calls the functions and reads the txt file
if __name__ == '__main__':

	print('E:/DO\ QTL\ MAPPING/')
	print("Input the above path if path needs to be specified since we need to use escape characters")

	# connect to database
	cnxn = createConnection('PARKSLAB', database)
	print('connected to the database: %s successfully!' %database)
	cursor = cnxn.cursor()

	# if create is specified in the command line
	if create:
		createTable(database, tablename)

	# check if path exists
	if not os.path.isdir(path):
   		print('path does not exist')
		exit(1)

    # change directory
	os.chdir(path)

	fileNames = []

	# for each of the files in the dir ending with txt, add to the list
	for file in os.listdir(path):
		if file.endswith(".gwas"):
			fileNames.append(file)

	for fileName in fileNames:

		print("Reading from:" + str(fileName))

		with open(fileName, 'r') as txtFile:

			txtReader = csv.reader(txtFile, delimiter = '\t')

			# Skip through the description header 
			fileFormat = next(txtReader)

			index = 0
			# Resort each row to resemble the database format
			for rows in txtReader:

				index = index + 1

				list = []

				list.append(fileName)
				list.append(rows[0])
				list.append(rows[3])
				list.append(rows[4])
				
				'''
				for i in range(0, len(list)):
					if list[i] == 'NA':
						list[i] = 'NULL'
				'''

				try:
						query = "insert into dbo.{!s}".format(tablename) +\
								"(trait, SNP, P, Beta)" + \
								" values ({!r}, {!r}, {!s}, {!s});".format(list[0], list[1], str(Decimal(list[2])), str(Decimal(list[3])))
						cursor.execute(query)
						cursor.commit()
						
				# write errmsg if file I/O exception
				except ValueError as ex:

					errmsg = "Warning: Value Error in " + str(fileName) + ", primary key is: {!r}".format(list[1]) + ", the line is: " + str(index)
					f = open("GC_{!r}_err.txt".format(tablename), "w")
					f.write(errmsg + "\n")
					f.write(str(list))
					f.close()

				except IndexError as iex :
					errmsg = "Warning: Index error in " + str(fileName) + ", primary key is: {!r}, {!r}".format(list[0], list[1])
					f = open("GC_{!r}_err.txt".format(tablename), "w")
					f.write(errmsg + "\n")
					f.write(str(list))
					f.close()

				except Exception as eex:
					print("Primary Key Integrity Violation for insert #" + str(index))

				else:
					print("Insert " + str(index) + " was successful!")
	
	cursor.commit()
	print("File Read Done!" + str(fileName))
	cnxn.close()
						