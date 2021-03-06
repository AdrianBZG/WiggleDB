Installing WiggleDB
===================

WiggleDB is a server application to run WiggleTools commands remotely.

Install WiggleTools
-------------------

```
git clone https://github.com/Ensembl/WiggleTools.git
cd WiggleTools
sudo sh ./easy_install.sh
chmod 755 bin/* python/wiggletools/*
cp -r python/wiggletools/ /usr/local/lib/python-2.7/dist-packages
cd ..
```
Install WiggleDB
----------------

```
git clone https://github.com/Ensembl/WiggleDB.git
cd WiggleDB
chmod 755 python/wiggleDB/* cgi/*
cp -r python/wiggletools/ /usr/local/lib/python-2.7/dist-packages
```

Prepare the data
----------------

1. Prepare a tab-delimited file, say datasets.tsv. The first five columns headers must be:	

	```
	location	name	type	annotation	assembly
	```
	Their content should be:
	- location: absolute path to the file
	- name: string (useful for annotations)
	- type: regions or signal
	- annotation: TRUE or FALSE
	- assembly: assembly name (e.g. GRCh38)
	- all other columns are your business. Try to select meaningful headers, and remove columns with too many values.

2. Prepare a tab-delimited file, say chromosome.lengths, with chromosome lengths and ensure it has 644 permissions, e.g.:

	```
	1	249250621
	10	135534747
	11	135006516
	12	133851895
	13	115169878
	14	107349540
	15	102531392
	```

3. Create an SQLite3 database:

	```	
	wiggleDB.py --database database.sqlite3 --load datasets.tsv
	wiggleDB.py --database database.sqlite3 --load assembly_name chromosome.lengths
	chmod 777 database.sqlite3
	```

4. Create a JSON file containing file attributes and allowed values:

	```
	wiggleDB.py --database database.sqlite3 --attributes > gui/datasets.attribs.json 
	chmod 644 gui/*
	```

5. Move the SQLite3 file to a location visible to all users.

Install AWS CLI
---------------

If you wish to push final results to S3 storage, create a config file with credentials, readable by all users.

Prepare the server
------------------

1. Create/choose log directory, ensure 777 permissions
2. Create/choose tmp directory, ensure 777 permissions
3. Create your own config file (see example in conf/wiggletools.conf) (Ensure 644 permissions, move to a visible location)
4. Copy the content of cgi/ to your Apache CGI directory, and edit the top of the CGI file, so that it points to your config file. 
	- Test by running wiggleCGI.py on the command line, without parameters
6. Copy the content of gui/ to your Apache web directory
	- check the URLs at the top of the Javascript file
