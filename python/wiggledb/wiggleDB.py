#!/usr/bin/env python 
# Copyright 2013 EMBL-EBI
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import re
import argparse
import sqlite3
import subprocess
import tempfile
import os
import os.path
import json
import math 
import numpy
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.patches as mpatches

verbose = False

###########################################
## Configuration file
###########################################

def read_config_file(filename):
	return dict(line.strip().split('\t') for line in open(filename) if line[0] != '#' and len(line) > 1)

###########################################
## Command line interface
###########################################

def normalise_spaces(string):
	if string is None:
		return None
	else:
		return re.sub("^\W+", "", re.sub("\W+", " ", string))

def get_options():
	parser = argparse.ArgumentParser(description='WiggleDB backend.')
	parser.add_argument('--db', '-d', dest='db', help='Database file',required=True)
	parser.add_argument('--wd', dest='working_directory', help='Data directory')
	parser.add_argument('-a',dest='a',help='A set of SQL constraints',nargs='*')
	parser.add_argument('-wa',dest='wa',help='WiggleTools command for A')
	parser.add_argument('-b',dest='b',help='A second set of SQL constraints',nargs='*')
	parser.add_argument('-wb',dest='wb',help='WiggleTools command for B')
	parser.add_argument('--wiggletools','-w',dest='fun_merge',help='Wiggletools command')
	parser.add_argument('--emails','-e',dest='emails',help='List of e-mail addresses for reminder',nargs='*')

	parser.add_argument('--load',dest='load',help='Datasets to load in database')
	parser.add_argument('--load_assembly',dest='load_assembly',help='Load assembly description with chromosome lengths', nargs=2)
	parser.add_argument('--load_annotations',dest='load_annots',help='References to load in database')
	parser.add_argument('--clean',dest='clean',help='Delete cached datasets older than X days', type=int)
	parser.add_argument('--cache',dest='cache',help='Dump cache info', action='store_true')
	parser.add_argument('--datasets',dest='datasets',help='Print dataset info', action='store_true')
	parser.add_argument('--user_datasets',dest='user_datasets',help='Print user dataset info', nargs='*')
	parser.add_argument('--clear_cache',dest='clear_cache',help='Reset cache info', nargs='*')
	parser.add_argument('--remember',dest='remember',help='Preserve dataset from garbage collection', action='store_true')
	parser.add_argument('--dry-run',dest='dry_run',help='Do not run the command, print wiggletools command', action='store_true')
	parser.add_argument('--result','-r',dest='result',help='Return status or end result of job', type=int)
	parser.add_argument('--attributes','-t',dest='attributes',help='Print JSON hash of attributes and values', action='store_true')
	parser.add_argument('--verbose','-v',dest='verbose',help='Turn on status output',action='store_true')
	parser.add_argument('--config','-c',dest='config',help='Configuration file')
	parser.add_argument('--annotations','-n',dest='annotations',help='Print list of annotation names', action='store_true')
	parser.add_argument('--jobs','-j',dest='jobs',help='Print list of jobs',nargs='*')
	parser.add_argument('--upload','-u',dest='upload',help='Upload dataset')
	parser.add_argument('--userid',dest='userid',help='User name')
	parser.add_argument('--description',dest='description',help='Uploaded dataset description',default='TEST')

	options = parser.parse_args()
	if all(X is None for X in [options.load, options.clean, options.result, options.datasets, options.clear_cache]) and not options.cache and not options.attributes and not options.annotations:
		assert options.a is not None, 'No dataset selection to run on'
		assert options.wa is not None, 'No dataset transformation to run on'
		if options.b is not None:
			assert options.fun_merge is not None, 'No action command (load,clean,compute) specified'

	options.wa = normalise_spaces(options.wa)	
	options.wb = normalise_spaces(options.wb)	
	options.fun_merge = normalise_spaces(options.fun_merge)

	if options.config is not None:
		config = read_config_file(options.config)
		if options.working_directory is None:
			options.working_directory = config['working_directory']
	else:
		config = None

	global verbose
	verbose = options.verbose
	return options, config

###########################################
## Convenience functions
###########################################

def run(cmd):
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	ret = p.wait()
	out, err = p.communicate()
	if ret != 0:
		if verbose:
			sys.stdout.write("Failed in running %s\n" % (cmd))
			sys.stdout.write("OUTPUT:\n%s\n" % (out))
			sys.stdout.write("ERROR:\n%s\n" % (out))
		raise BaseException
	return out

###########################################
## Creating a database
###########################################

def create_database(cursor, filename):
	if verbose:
		print 'Creating database'
	create_cache(cursor)
	create_dataset_table(cursor, filename)
	create_annotation_dataset_table(cursor)
	create_user_dataset_table(cursor)
	create_assembly_table(cursor)
	create_chromosome_table(cursor)

def create_cache(cursor):
	cursor.execute('''
	CREATE TABLE IF NOT EXISTS
	cache
	(
	merge varchar(100),
	mergeA varchar(100),
	filesA varchar(10000),
	mergeB varchar(100),
	filesB varchar(10000),
	location varchar(1000) UNIQUE,
	userid varchar(100),
	last_query datetime
	)
	''')

def create_dataset_table(cursor, filename):
	file = open(filename)
	items = file.readline().strip().split('\t')
	assert items[0] == 'location', "Badly formed dataset table, please ensure the first column refers to location:\n" + items[0]
	assert items[1] == 'id', "Badly formed dataset table, please ensure the first column refers to location:\n" + items[0]
	header = '''
			CREATE TABLE IF NOT EXISTS 
			datasets 
			(
			location varchar(1000),
			id string(100),
		 '''
	cursor.execute('\n'.join([header] + [",\n".join(['%s varchar(255)' % X for X in items[2:]])] + [')']))

	cursor.execute('SELECT * FROM datasets').fetchall()
	column_names = [X[0] for X in cursor.description]
	assert column_names == items, 'Mismatch between the expected columns: \n%s\nAnd the columns in file:\n%s' % ("\t".join(column_names), '\t'.join(items))

	for line in file:
		cursor.execute('INSERT INTO datasets VALUES (%s)' % ",".join("'%s'" % X for X in line.strip().split('\t')))
	file.close()

def create_annotation_dataset_table(cursor):
	header = '''
			CREATE TABLE IF NOT EXISTS 
			annotation_datasets 
			(
			name varchar(1000) UNIQUE,
			location varchar(1000),
 			description varchar(100),
			count integer
			)
		 '''
	cursor.execute(header)

def create_user_dataset_table(cursor):
	header = '''
			CREATE TABLE IF NOT EXISTS 
			user_datasets 
			(
 			name varchar(100),
			location varchar(1000),
			userid varchar(100),
			count integer,
			has_history integer
			)
		 '''
	cursor.execute(header)

def create_chromosome_table(cursor):
	header = '''
			CREATE TABLE IF NOT EXISTS 
			chromosomes
			(
 			name varchar(1000) UNIQUE
			)
		 '''
	cursor.execute(header)

def create_assembly_table(cursor):
	header = '''
			CREATE TABLE IF NOT EXISTS 
			assemblies
			(
 			name varchar(100) UNIQUE,
 			location varchar(1000)
			)
		 '''
	cursor.execute(header)

###########################################
## Garbage cleaning 
###########################################

def clean_database(cursor, days):
	for location in cursor.execute('SELECT location FROM cache WHERE julianday(\'now\') - julianday(last_query) > %i' % days).fetchall():
		if verbose:
			print 'Removing %s' % location[0]
		if os.path.exists(location[0]):
			os.remove(location[0])
	cursor.execute('DELETE FROM cache WHERE julianday(\'now\') - julianday(last_query) > %i' % days)

###########################################
## Dataset attributes
###########################################

def get_dataset_attributes_2(cursor):
	return [X[1] for X in cursor.execute('PRAGMA table_info(datasets)').fetchall()]

def get_dataset_attributes(cursor):
	return list(set(get_dataset_attributes_2(cursor)) - set(["location","id"]))

def get_attribute_values_2(cursor, attribute):
	return [X[0] for X in cursor.execute('SELECT DISTINCT %s FROM datasets' % (attribute)).fetchall()]

def get_attribute_values(cursor):
	return dict((attribute, get_attribute_values_2(cursor, attribute)) for attribute in get_dataset_attributes(cursor))

def get_annotations(cursor):
	return cursor.execute('SELECT * FROM annotation_datasets').fetchall()

def get_datasets(cursor):
	res = cursor.execute('SELECT * FROM datasets').fetchall()
	return [[X[0] for X in cursor.description]] + res

def attribute_selector(attribute, params):
	return "( %s )" % " OR ".join("%s=:%s_%i" % (attribute,attribute,index) for index in range(len(params[attribute])))

def denormalize_params(params):
	return dict(("%s_%i" % (attribute, index),value) for attribute in params for (index, value) in enumerate(params[attribute]))

def get_dataset_locations(cursor, params):
	# Quick check that all the keys are purely alphanumeric to avoid MySQL injections
	assert not any(re.match('\W', X) is not None for X in params)
	query = " AND ".join(attribute_selector(X, params) for X in params)
	if len(query) > 0:
		query = ' WHERE ' + query 
	if verbose:
		print 'Query: SELECT location FROM datasets' + query
		print 'Where:' + str(denormalize_params(params))
	res = cursor.execute('SELECT location FROM datasets' + query, denormalize_params(params)).fetchall()
	if verbose:
		print 'Found:\n' + "\n".join(X[0] for X in res)
	return sorted(X[0] for X in res)

###########################################
## Annotation files
###########################################

def get_annotation_dataset_locations(cursor, names, userid = None):
	locations = []
	for name in names:
		res = cursor.execute('SELECT location FROM annotation_datasets WHERE name=?', (name,)).fetchall()
		if len(res) > 0:
			locations.append(res[0][0])
		elif userid is not None:
			res = cursor.execute('SELECT location FROM user_datasets WHERE name=? AND userid=?', (name,userid)).fetchall()
			if len(res) > 0:
				locations.append(res[0][0])
	return locations

def get_annotation_dataset_description(cursor, name):
	res = cursor.execute('SELECT description FROM annotation_datasets WHERE name=?', (name,)).fetchall()
	if len(res) > 0:
		return {'name':name, 'description':res[0][0]}
	else:
		return {'ERROR'}

###########################################
## Assemblies and chromosomes
###########################################

def register_new_chromosomes(cursor, location):
	file = open(location)
	for line in file:
		items = line.strip().split('\t')
		if len(items) == 2:
			cursor.execute('INSERT INTO chromosomes (name) VALUES (?)', (items[0], ))
	file.close()

def register_assembly(cursor, name, location):
	cursor.execute('INSERT INTO assemblies (name, location) VALUES (?,?)', (name, location))
	register_new_chromosomes(cursor, location)

def get_assembly_file(cursor):
	res = cursor.execute('SELECT location FROM assemblies').fetchall()
	if len(res) > 0:
		return res[0][0]
	else:
		raise
def get_chromosomes(cursor):
	return [X[0] for X in cursor.execute('SELECT name FROM chromosomes').fetchall()]

###########################################
## Annotations
###########################################

def count_regions(location):
	if location[-4:] == '.bed' or location[-4:] == '.txt':
		return int(run('cat %s | wc -l' % location).strip())
	elif location[-3:] == '.bb':
		return int(run('bigBedToBed %s stdout | grep -v ^# | wc -l' % (location)).strip())
	else:
		raise BaseException

def add_annotation_dataset(cursor, location, name, description):
	cursor.execute('INSERT INTO annotation_datasets (name,location,description, count) VALUES (?,?,?,?)', (name, location, description, count_regions(location)))

def get_user_dataset_locations(cursor, names, userid):
	locations = []
	for name in names:
		res = cursor.execute('SELECT location FROM user_datasets WHERE name=? AND userid=?', (name,userid)).fetchall()
		if len(res) > 0:
			locations.append(res[0][0])
		else:
			return []
	return locations

###########################################
## General search 
###########################################

def get_locations(cursor, params, userid):
	if 'annot_name' in params:
		return get_annotation_dataset_locations(cursor, params['annot_name'], userid)
	elif 'user_name' in params:
		return get_user_dataset_locations(cursor, params['user_name'], userid)
	else:
		return get_dataset_locations(cursor, params)

###########################################
## Search cache
###########################################

def reset_time_stamp(cursor, location):
	cursor.execute('UPDATE cache SET last_query= date(\'now\') WHERE location = \'%s\'' % location)

def get_precomputed_location(cursor, merge, mergeA, filesA, mergeB, filesB):
	filesA_str = " ".join(filesA)
	filesB_str = " ".join(filesB)
	reports = cursor.execute('SELECT location FROM cache WHERE merge = ? AND mergeA = ? AND filesA = ? AND mergeB = ? AND filesB = ?', (merge, mergeA, filesA_str, mergeB, filesB_str)).fetchall()
	if len(reports) > 0:
		if verbose:
			print 'Found pre-computed file for query: %s' % cmd
			print reports[0]
		location = reports[0][0]
		reset_time_stamp(cursor, location)
		return location
	else:
		if verbose:
			print 'Did not find pre-computed file for query: %s' % cmd
		return None

def is_tracked(cursor, location):
	return cursor.execute('SELECT COUNT(*) FROM cache WHERE location = ?', (location,)).fetchall()[0][0] > 0

def copy_to_longterm(data, config):
	if 's3_bucket' in config:
		os.environ['AWS_CONFIG_FILE'] = config['aws_config']
		run("aws s3 cp %s s3://%s/%s --acl public-read" % (data, config['s3_bucket'], os.path.basename(data)))

def get_annotation_counts(cursor, userid):
	counts = [(X[0], X[1]) for X in cursor.execute('SELECT name, count FROM annotation_datasets').fetchall()]
	if userid is not None:
		counts += [(X[0], X[1]) for X in cursor.execute('SELECT name, count FROM user_datasets WHERE userid=?', (userid,)).fetchall()]
	return dict(counts)

def make_barchart(counts, total, rev_counts, totals, labels, out, format='pdf'):
	width = 0.35

	ind = numpy.arange(len(labels)) - width/2 
	heights = [X/float(total) for X in counts]
	assert all(X <= 1 and X >= 0 for X in heights), (counts, total, heights)
	errors = [ math.sqrt(2*X*(1-X) / total) for X in heights ]
	pyplot.bar(ind, heights, width, yerr=errors, color='b')

	ind = numpy.arange(len(labels)) + width/2 
	heights = [X/float(Y) for X, Y in zip(rev_counts, totals)]
	assert all(X <= 1 and X >= 0 for X in heights), (zip(rev_counts, totals), heights)
	errors = [ math.sqrt(2*X*(1-X) / Y) for X, Y in zip(heights, totals) ]
	pyplot.bar(ind, heights, width, yerr=errors, color='r')

	blue_patch = mpatches.Patch(color='blue', label='Precision')
	red_patch = mpatches.Patch(color='red', label='Recall')
	pyplot.legend((blue_patch, red_patch), ('Specificity', 'Sensitivity'))

	pyplot.xticks(ind, labels)
	pyplot.savefig(out, format=format)

def launch_quick_compute(conn, cursor, fun_merge, fun_A, data_A, fun_B, data_B, options, config):
	cmd_A = " ".join([fun_A] + data_A + [':'])

	if len(data_B) > 0:
		merge_words = fun_merge.split(' ')

		assert fun_merge is not None
		if fun_B != "":
			cmd_B = " ".join([fun_B] + data_B + [':'])
		else:
			cmd_B = " ".join(data_B)

		if merge_words[0] == 'overlaps':
			assert "annot_name" in options.b
			fh, destination = tempfile.mkstemp(suffix='.txt',dir=options.working_directory)
			total = int(run("wiggletools write_bg - %s | wc -l" % cmd_A))
			if total > 0:
				counts = []
				for annotation in data_B:
					counts.append(int(run('wiggletools write_bg - overlaps %s %s | wc -l' % (annotation, cmd_A)).strip()))
			
				rev_counts = []
				for annotation in data_B:
					rev_counts.append(int(run('wiggletools write_bg - overlaps %s %s | wc -l' % (cmd_A, annotation)).strip()))

				annotation_count_dict = get_annotation_counts(cursor, options.userid)
				annotation_counts = []
				for annotation in options.b['annot_name']:
					annotation_counts.append(annotation_count_dict[annotation])

				out = open(destination, "w")
				for name, count, annotation_count in zip(options.b['annot_name'], counts, annotation_counts):
					out.write("\t".join(map(str, [name, count, annotation_count])) + "\n")
				out.write("\t".join(['ALL', str(total)]) + "\n")
				out.close()

				make_barchart(counts, total, rev_counts, annotation_counts, options.b['annot_name'], destination + '.png', format='png')
		else:
			fh, destination = tempfile.mkstemp(suffix='.bed',dir=options.working_directory)
			os.remove(destination)
			run(" ".join(['wiggletools','write_bg', destination, fun_merge, cmd_A, cmd_B]))
	else:
		fh, destination = tempfile.mkstemp(suffix='.bed',dir=options.working_directory)
		os.remove(destination)
		run(" ".join(['wiggletools','write_bg', destination, cmd_A]))

	return destination

def make_normalised_form(fun_merge, fun_A, data_A, fun_B, data_B):
	cmd_A = funA + ":" + " ".join(data_A)
	if data_B is not None:
		if fun_B is not None:
			cmd_B = funB + ":" + " ".join(data_B)
		else:
			cmd_B = " ".join(data_B)
		res = ";".join([fun_merge, cmd_A, cmd_B])
	else:
		res = cmd_A

	if verbose:
		print 'CMD A: ' + cmd_A
		print 'CMD B: ' + str(cmd_B)
		print 'CMD  : ' + res

	return res

def form_filter(cursor, f, userid):
	if len(f) != 3:
		return ""
	elif f[1] == 0:
		return f[0] + " " + get_annotation_dataset_locations(cursor, [f[2]], userid)[0]
	else:
		return " ".join([f[0],"extend",f[1],get_annotation_dataset_locations(cursor, [f[2]], userid)[0]])

def form_filters(cursor, filters, userid):
	if filters is not None:
		return " ".join([form_filter(cursor, X, userid) for X in filters])
	else:
		return ""

def request_compute(conn, cursor, options, config):
	if options.fun_merge is None:
		options.fun_merge = ""
	fun_A = form_filters(cursor, options.filters_a, options.userid) + " " + options.wa
	data_A = get_locations(cursor, options.a, options.userid)
	if len(data_A) == 0:
		 return {'status':'INVALID'}

	if options.b is not None:
		if options.wb is not None:
			fun_B = form_filters(cursor, options.filters_b, options.userid) + " " + options.wb
		else:
			fun_B = "" 
		data_B = get_locations(cursor, options.b, options.userid)
		options.countB = len(data_B)
		if len(data_B) == 0:
			 return {'status':'INVALID'}
	else:
		data_B = []
		fun_B = "" 

	prior_result = get_precomputed_location(cursor, options.fun_merge, fun_A, data_A, fun_B, data_B)
	if prior_result is None:
		destination = launch_quick_compute(conn, cursor, options.fun_merge, fun_A, data_A, fun_B, data_B, options, config)
		if os.stat(destination).st_size == 0:
			return {'status':'EMPTY'}
		else:
			copy_to_longterm(destination, config)
			if destination[-4:] == ".txt":
				copy_to_longterm(destination + ".png", config)
			if destination[-4:] == ".bed":
				run(" ".join(['wigToBigWig',destination, get_assembly_file(cursor), destination + '.bw']))
		filesA_str = " ".join(data_A)
		if data_B is None:
			filesB_str = None
		else:
			filesB_str = " ".join(data_B)
		cursor.execute('INSERT INTO cache (merge,mergeA,filesA,mergeB,filesB,location,userid,last_query) VALUES (?,?,?,?,?,?,?,date("now"))', (options.fun_merge, fun_A, filesA_str, fun_B, filesB_str, destination, options.userid))
		return {'location': destination, 'status':'DONE'}
	else:
		reset_time_stamp(cursor, prior_result)
		return {'location': prior_result, 'status':'DONE'}

###########################################
## Upload file
###########################################

def fetch_user_datasets(cursor, userid):
	return cursor.execute('SELECT * FROM user_datasets WHERE userid=?', (userid,)).fetchall()

def get_user_datasets(cursor, userid):
	return {'files': dict((X[0], X[1]) for X in cursor.execute('SELECT name, has_history FROM user_datasets WHERE userid=?', (userid,)).fetchall())}

def remove_user_datasets(cursor, name, userid):
	locations = get_user_dataset_locations(cursor, [name], userid)
	cursor.execute('DELETE FROM user_datasets WHERE name=? AND userid=?', (name, userid))
	return {'status':'SUCCESS'}

def wget_dataset(url, dir):
	fh, destination = tempfile.mkstemp(suffix="." + url.split('.')[-1],dir=dir)
	run("wget %s -O %s" % (url, destination))
	return destination

def check_file_integrity(file, cursor):
	suffix = file.split('.')[-1]
	if suffix == 'bed' or suffix == 'txt':
		try:
			chromosomes = get_chromosomes(cursor)
			bad_assembly = True
			found_chromosomes = set()
			elems = []
			fh = open(file, "r")
			for line in fh:
				if len(line.strip()) == 0 or line[0] == '#' or line[0] == '-':
					continue
				items = line.strip().split('\t')
				if len(items) < 3:
					return {'status':'MALFORMED_INPUT', 'format':'Bed', 'line': line}
				if items[0] in chromosomes:
					bad_assembly = False
				found_chromosomes.add(items[0])
				elems.append((items[0], int(items[1]), int(items[2])))
			fh.close()
			if bad_assembly:
				return {'status':'WRONG_ASSEMBLY', 'found_chromosomes':list(found_chromosomes), 'expected_chromosomes': chromosomes}
			fh = open(file, "w")
			for elem in sorted(elems):
				fh.write("\t".join(map(str, elem)) + "\n")
			fh.close()
			return
		except:
			return {'status':'MALFORMED_INPUT', 'format':'Bed'}
	elif suffix == 'bb':
		try:
			run("bigBedInfo %s" % (file))
		except:
			return {'status':'MALFORMED_INPUT', 'format':'bigBed'}
		
	return {'status':'MALFORMED_INPUT','format':'unrecognized'}

def register_user_dataset(cursor, file, description, userid, has_history=False): 
	cursor.execute('INSERT INTO user_datasets (name,location,userid,count,has_history) VALUES (?,?,?,?,?)', (description, file, userid, count_regions(file), int(has_history)))

def name_already_used(cursor, description, userid):
	return len(get_user_dataset_locations(cursor, [description], userid)) > 0 or len(get_annotation_dataset_locations(cursor, [description])) > 0

def upload_dataset(cursor, dir, url, description, userid):
	if name_already_used(cursor, description, userid):
		return {'status':'NAME_USED','name':description}

	try:
		file = wget_dataset(url, dir);
		if file is None:
			raise
		ret = check_file_integrity(file, cursor)
		if ret is not None:
			return ret
		register_user_dataset(cursor, file, description, userid)
		return {'status':'UPLOADED','name':description}
	except:
		raise
		return {'status':'UPLOAD_FAILED','url':url}

def save_dataset(cursor, dir, fileitem, description, userid):
	if name_already_used(cursor, description, userid):
		return {'status':'NAME_USED','name':description}

	try:
		fh, destination = tempfile.mkstemp(suffix=".bed",dir=dir)
		file = open(destination, "w")
		file.write(fileitem.file.read())
		file.close()
		ret = check_file_integrity(destination, cursor)
		if ret is not None:
			return ret
		register_user_dataset(cursor, destination, description, userid)
		return {'status':'UPLOADED','name':description}
	except:
		raise
		return {'status':'UPLOAD_FAILED','url':url}

def reassign_dataset(cursor, filename, description, userid):
	if name_already_used(cursor, description, userid):
		return {'status':'NAME_USED','name':description}
	else:
		register_user_dataset(cursor, filename, description, userid, has_history=True)
		return {'status':'UPLOADED','name':description, 'has_history':1}

###########################################
## History
###########################################

def get_user_annotation_history(cursor, location, userid):
	res = cursor.execute('SELECT name FROM user_datasets WHERE location=? AND userid = ?', (location,userid)).fetchall()
	if len(res) > 0:
		return "<LOCAL_DATA:%s>" % res[0][0]
	else:
		return None

def get_dataset_history(cursor, location):
	res = cursor.execute('SELECT id FROM datasets WHERE location=?', (location,)).fetchall()
	if len(res) > 0:
		return "<LOCAL_DATA:%s>" % res[0][0]
	else:
		return None
	
def get_annotation_history(cursor, location):
	res = cursor.execute('SELECT merge, mergeA, filesA, mergeB, filesB, userid FROM cache WHERE location=?', (location,)).fetchall()
	if len(res) > 0:
		merge, mergeA, filesA, mergeB, filesB, userid = res[0]
		if filesA != "":
			filesA = " ".join([get_location_history(cursor, X, userid) for X in filesA.split(" ")]) + " :"
		if filesB != "":
			filesB = " ".join([get_location_history(cursor, X, userid) for X in filesB.split(" ")]) + " :"
		return " ".join([merge, mergeA, filesA, mergeB, filesB])
	else:
		return None

def get_location_history(cursor, location, userid):
	dataset_string = get_dataset_history(cursor, location)
	if dataset_string is not None:
		return dataset_string

	annotation_string = get_annotation_history(cursor, location)
	if annotation_string is not None:
		return annotation_string

	user_string = get_user_annotation_history(cursor, location, userid)
	if user_string is not None:
		return user_string

	return "<UNKNOWN:%s>" % location

def get_named_file_history(cursor, name, userid):
	res = cursor.execute('SELECT location, has_history, userid FROM user_datasets WHERE name=? AND userid=?', (name,userid)).fetchall()
	if len(res) > 0:
		location = res[0][0]
		has_history = res[0][1]
		if has_history:
			return get_location_history(cursor, location, userid)
		else:
			return "<USER_DATA:%s>" % name
	else:
		return None

###########################################
## Main
###########################################

def main():
	options, config = get_options()
	conn = sqlite3.connect(options.db)
	cursor = conn.cursor()

	if options.load is not None:
		create_database(cursor, options.load)
	elif options.load_annots is not None:
		for line in open(options.load_annots):
			items = line.strip().split('\t')
			assert len(items) == 3
			add_annotation_dataset(cursor, items[0], items[1], items[2])
	elif options.load_assembly is not None:
		name, location = options.load_assembly
		register_assembly(cursor, name, os.path.abspath(location))
	elif options.clean is not None:
		clean_database(cursor, options.clean)
	elif options.cache:
		for entry in cursor.execute('SELECT * FROM cache').fetchall():
			print entry
	elif options.clear_cache is not None:
		if len(options.clear_cache) == 0:
			cursor.execute('DROP TABLE cache')
			create_cache(cursor)
			cursor.execute('DROP TABLE user_datasets')
			create_user_dataset_table(cursor)
		else:
			remove_jobs(cursor, options.clear_cache)
	elif options.attributes:
		print json.dumps(get_attribute_values(cursor))
	elif options.datasets:
		print "\n".join("\t".join(map(str, X)) for X in get_datasets(cursor))
	elif options.annotations:
		print "\n".join("\t".join(map(str, X)) for X in get_annotations(cursor))
	elif options.upload is not None:
		print json.dumps(upload_dataset(cursor, options.working_directory, options.upload, options.description, options.userid))
	elif options.user_datasets is not None:
		if len(options.user_datasets) == 0:
			print "\n".join("\t".join(map(str, X)) for X in cursor.execute("SELECT * FROM user_datasets").fetchall())
		else:
			for user in options.user_datasets:
				print json.dumps(fetch_user_datasets(cursor, user))
	else:
		if options.a is not None:
			res = dict()
			for constraint in options.a:
				attribute, value = constraint.split("=")
				if attribute not in res:
					res[attribute] = []
				res[attribute].append(value)
			options.a = res
		if options.b is not None:
			res = dict()
			for constraint in options.b:
				attribute, value = constraint.split("=")
				if attribute not in res:
					res[attribute] = []
				res[attribute].append(value)
			options.b = res
		print json.dumps(request_compute(conn, cursor, options, config))

	conn.commit()
	conn.close()

if __name__=='__main__':
	main()
