import pymysql
import configparser
import logging
import os
import re
import warnings
import sys
import zipfile
from datetime import datetime
import glob
import subprocess

from hashlib import md5

def get_md5_of_zip(file_name):

    m = md5()
    with open(file_name, "rb") as f:
        data = f.read() #read file in chunk and call update on each chunk if file is large.
        m.update(data)
        return m.hexdigest()

archives = ['archive/a2018/FASTQ']
scratch_base = '/nfs/fastq_temp2/CLEAN_ME/'
#scratch_ext = ['EXOME','EXOME_1.7','EXOME_ssd','GENOME','GENOME_1.7','GENOME_ssd','RNASEQ']
scratch_ext = ['EXOME']

def get_connection_mml(database):
    try:
        reader = configparser.RawConfigParser()
        reader.read('/home/mml2204/.my.cnf')
        db = 'client' + database
        db_host = reader.get(db, 'host')
        db_user = reader.get(db, 'user')
        db_pass = reader.get(db,'password')
        connection = pymysql.connect(host=db_host,user=db_user,passwd=db_pass,db='sequenceDB')
        return connection
    except pymysql.Error:
        traceback.print_exc()
        sys.exit("Wrong username/database or password found, please try again")

def run_query_mml(query,connection):
    try:
            cursor = connection.cursor()
            cursor.execute(query)
            results = cursor.fetchall()
            connection.commit()
            return results, cursor.rowcount
    except pymysql.Error:
        print("Problem with database connection. Query %s didnt execute" %query)
        sys.exit("Error")

def subprocess_exec(cmd):
    p1 = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = p1.communicate()
    if err.decode() != '':
        print("problem with {} out:{}, err:{}".format(cmd,out.decode().strip(),err.decode().strip()))
    return out.decode().strip().split('\n')

def main():
  for item in scratch_ext:
    if 'EXOME' in item:
        seqtype = 'EXOME'
    elif 'GENOME' in item:
        seqtype = 'GENOME'
    elif 'RNASEQ' in item:
        seqtype = 'RNASEQ'
    else:
        raise Exception("invalid seqtype")
    ignore = []
    if len(sys.argv) == 2:
        with open(sys.argv[1],'r') as f:
            for line in f:
                ignore.append(line.strip())

    full_path = scratch_base + item
    samps = os.listdir(full_path)
    connection = get_connection_mml('sequenceDB')
    for sample in samps:
        cmd = "SELECT status from prepT where chgvid='{}' and sample_type='{}'".format(sample,seqtype)
        info,num_rows = run_query_mml(cmd,connection)
        if num_rows > 1 and ( info[0][0] == 'In DragenDB' or info[1][0] == 'In DragenDB' ):
            pass
        elif num_rows != 1:
            print("TOO MANY ROWS RETURNED")
            print(cmd)
            print(info,num_rows)
            continue
        if num_rows == 1 and info[0][0] != 'In DragenDB':
            print("not in DragenDB {}, {}".format(sample,seqtype))
            continue
        else:
            if '{}/{}'.format(full_path,sample) in ignore:
                print("IGNORE THIS ONE {}/{}".format(full_path,sample))
                continue
            flos = os.listdir(full_path+'/'+sample+'/')
            for flowcell in flos:
                #find flowcell in archive location
                print('{}/{}/{}'.format(full_path,sample,flowcell))
                if '{}/{}/{}'.format(full_path,sample,flowcell) in ignore:
                    print("IGNORE THIS ONE {}/{}/{}".format(full_path,sample,flowcell))
                    continue
                for archive in archives:
                   archive_loc = '/nfs/{}/{}/{}/{}'.format(archive,seqtype,sample,flowcell)
                   temp_path = '{}/{}/{}'.format(full_path,sample,flowcell)
                   #print(archive_loc)
                   #print(temp_path)
                   if os.path.exists(archive_loc):
                       archive_files = os.listdir(archive_loc)
                       temp_files = os.listdir(temp_path)
                       for filename in temp_files:
                           if os.path.islink('{}/{}'.format(temp_path,filename)):
                               print('{}/{} is a symlink'.format(temp_path,filename))
                               continue
                           if filename not in archive_files:
                               print('{} does not exist in archive {}/{}'.format(filename,archive_loc,filename))
                               continue
                           elif os.path.isdir('{}/{}'.format(temp_path,filename)):
                               print('{}/{} not a file.IGNORING'.format(temp_path,filename))
                               continue
                           else:
                               t_size = os.path.getsize('{}/{}'.format(temp_path,filename))
                               ar_size = os.path.getsize('{}/{}'.format(archive_loc,filename))
                               print("{}/{}".format(temp_path,filename)) 
                               #temp_md5 = get_md5_of_zip("{}/{}".format(temp_path,filename))
                               #ar_md5 = get_md5_of_zip("{}/{}".format(archive_loc,filename))
                               #print(temp_md5)
                               #print(ar_md5)
                               #if  t_size!= ar_size or temp_md5 != ar_md5:
                               if  t_size!= ar_size:
                                   warnings.warn("size not equal: ({0},{1}/{2}), ({3},{4}/{2})".format(t_size,temp_path,filename,ar_size,archive_loc))
                                   #warnings.warn("checksums not equal: ({0},{1}/{2}), ({3},{4}/{2}) md5sums not equal: ({5},{1}/{2}), ({6},{4}/{2})".format(t_size,temp_path,filename,ar_size,archive_loc, temp_md5, ar_md5))
                                   continue
                               else:
                                   cmd = "rsync -zpvv --remove-source-files --no-owner --no-perms {0}/{1} {2}/{1}".format(temp_path,filename,archive_loc)  
                                   #cmd="rm -v {}/{}".format(temp_path,filename)
                                   print(cmd)
                                   subprocess_exec(cmd)
                    
main()
