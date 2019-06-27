from utils import *
import sys
import argparse
import csv
import os
import pprint

"""
INPUT:
1) the "manifest" file with an optional additional column for 'file_name' which contains pipe separated list of filenames
    #if not manifest file is provided, sift through all 'External Data...' samples in DB
2) optional directory location. if not given it checks:
    a) /nfs/tx/tx_*
    b) /nfs/seqscratch10/tx_temp/tx_*
    c) /nfs/seqscratch09/tx_temp/tx_*
    d) /nfs/seqscratch11/tx_temp/tx_*
    e) /nfs/seqscratch12/tx_temp/tx_*

#also add in a check to test if bam / fastq / cram is not truncated.
#also add archival to /nfs/archive/p2018/FASTQ/EXOME/<sample_internal_name/<fcillumid>/<fastq.gz> or /nfs/archive/p2018/FASTQ/EXTERNAL_BAMS/EXOME/<sample_internal_name>/<bam or cram file>
#verify sample_type = 'EXOME'
#if exome, exomekit != 'N/A'
#also do some file sz checks

OUTPUT:
1)mapping from  chgvids to fastq, bam, cram files in directory; if cram mapping must contain a reference file.
2)unmapped fastq, bam, cram files - 2 categories 'map' and 'possible_map'.
3)unmapped  chgvids
"""

file_types = ['fastq','bam','cram']
ref_genome_extensions = ('.fa','.fai','.fa.gz','.fa.bgz')


default_file_dir = ['/nfs/tx/']#can add more later 
#imp_fields_from_man = ['CHGVID','OrigID','AKA','SeqType','SubProject','Protocol','CurrProjLeader','GAFbin','FamilyID','FamilyRelationProband','exomeKit','FundCode']
imp_fields_from_man = ['CHGVID','OrigID','AKA','SeqType','exomeKit']
def get_all_files_in_dir(file_basedirpath):
    
    all_files = dict((file_type,[]) for file_type in file_types)#note how fromkeys() doesnt work: https://stackoverflow.com/questions/15516413/dict-fromkeys-all-point-to-same-list
    all_files['ref'] = []
    all_files['other'] = []

    with open(file_basedirpath,'r') as basedir_files:
        for line in basedir_files:
            line = line.strip()
            if not os.path.exists(line):
                continue
            basedirpath = line
            for dirpath, dirnames, filenames in os.walk(basedirpath):
                for f in filenames:
                    found = 0 
                    if f.endswith(ref_genome_extensions):
                        all_files['ref'].append(os.path.join(dirpath,f))
                        found = 1 
                    else:    
                        for ext in file_types:
                            if f.endswith('.' + ext) or f.endswith('.' + ext + '.gz'):
                                full_filename = os.path.join(dirpath,f)
                                all_files[ext].append(full_filename)
                                found = 1 
                                break
                    if found == 0:
                        all_files['other'].append(os.path.join(dirpath,f))
    
    pprint.pprint(all_files)
    return(all_files)


def find_substr(base_str,sub_str):
    #print(base_str,sub_str)
    sub_pos = 0
    found_pos = []
    start = -1
    for base_pos in range(len(base_str)):
        #print(base_pos,base_str[base_pos],sub_str[sub_pos])
        if base_str[base_pos] in ['.','_','-']:
            if sub_pos == len(sub_str):
                #print("here 1")
                return 1
            if sub_str[sub_pos] in ['.','_','-']:
                base_pos += 1
                sub_pos += 1
                #print("here 2")
                continue
            #if sub_pos == len(sub_str):
            #    return 1
        elif sub_pos == len(sub_str):
            #print("here 3")
            return 0
        if base_str[base_pos].lower() == sub_str[sub_pos].lower():#for now keep it case sensitive, consider making case insensitive later ? 
            #if sub_pos == 0:
            #    start = base_pos
            #print("here 4")
            sub_pos += 1
        else:
            #print("here 5")
            sub_pos = 0
        base_pos += 1


def add_to_map(chgvid,origid,sure_mapping,prob_mapping,mapped_files,file_type,file_dict):
        for filename in file_dict[file_type]:
            #Ambry file structure: 17_461959_S4_L001_R2_001.fastq.gz and the origid would be 17-461959
            #for matches to 17_461959.bam
            #matches to origd
            #for origid[count] in file_to_map
    
            #add support for the UDN's later and also bam / cram support
    
            filename_split = filename.split('/')
            if find_substr(filename_split[-1],origid) or find_substr(filename_split[-1],chgvid):
                mapped_files[file_type].append(filename)
                if chgvid not in sure_mapping:
                    sure_mapping[chgvid] = {file_type:[filename]}
                else:
                    if file_type not in sure_mapping[chgvid]:
                        sure_mapping[chgvid][file_type] = [filename]
                    else:
                        sure_mapping[chgvid][file_type].append(filename)
                #add sanity checks for files
            elif filename_split[-1].find(chgvid) == 1 or filename_split[-1].find(origid) == 1:
                mapped_files[file_type].append(filename)
                if chgvid not in prob_mapping:
                    prob_mapping[chgvid] = {file_type:[filename]}
                else:
                    if file_type not in prob_mapping[chgvid]:
                        prob_mapping[chgvid][file_type] = [filename]
                    else:
                        prob_mapping[chgvid][file_type].append(filename)    

def get_file_data_mapping(chg,aka,orig,dup_chg,dup_orig,dup_aka,file_dict):
    sure_mapping = {}
    prob_mapping = {}
    dup_mapping = {} 
    unmapped = {}
    unmapped_chg = []
    #in_db = dict((file_type,[]) for file_type in file_types)
    mapped_files = dict((file_type,[]) for file_type in file_types)

    for count in range(len(chg)):
        print("sample number: {}".format(count))
        print(chg[count])
        print(orig[count])
        #print(dup_orig)
        if chg[count] in dup_chg or ( orig[count] != 'N/A' and orig[count] != '' and orig[count] in dup_orig ):
                dup_mapping[chg[count]] = "chg/orig_duplicated_in_manifest"
                continue
        for file_type in file_types:
            add_to_map(chg[count],orig[count],sure_mapping,prob_mapping,mapped_files,file_type,file_dict)
            
            #if using aka would need to consider ';' add check for aka later
            
            #not found in filename, check dir
    for file_type in mapped_files:
        unmapped[file_type] = list(set(file_dict[file_type]) - set(mapped_files[file_type]))
    
    """
    #if unampped are already in DragenDB, remove from unmapped and add to in DragenDB dict
    unmapped_list = []
    for file_type in unmapped:
        for filename in unmapped[file_type]:
            unmapped_list.extend(filename.split('/'))
    unmapped_list = ["'" + entry + "'" for entry in set(unmapped_list)]
    conn_sdb = get_connection_mml('sequenceDB','sequenceDB')
    #cmd = "select sample_internal_name as chgvid,status  from prepT where sample_internal_name in {}".format(','.join(unmapped_list))
    unmapped_found, numrows = run_query_mml_no_commit(cmd,conn_sdb)
    conn_sdb.close()
    for entry in unmapped_found:
        for file_type in unmapped:
            for filename in unmapped[file_type]:
                if '/' + entry['chgvid'] + '/' in filename:
    """
    
    
    for entry in chg:
        if entry not in sure_mapping.keys() and  entry not in prob_mapping.keys():
            unmapped_chg.append(entry)

    for sample in sure_mapping:
        if 'fastq' in sure_mapping[sample]:
            if len(sure_mapping[sample]['fastq']) %2 != 0:
                print("!!!WARNING!!!: missing data files for {}. Following files are present {}".format(sample,'\n'.join(sure_mapping[sample]['fastq'])))
                       
    print("sure")
    pprint.pprint(sure_mapping)
    print("prob")
    print(prob_mapping)
    print("dup")
    print(dup_mapping)
    print("unmapped")
    pprint.pprint(unmapped)
    print("mapped")
    pprint.pprint(mapped_files)
    print("unmapped chgvid")
    pprint.pprint(unmapped_chg)



def get_recent_ext_sub(conn_sdb):    
    #get all samples with prepT.status = 'External Data Submitted', check their chgvids, origids, akas for match with "unmapped" files in this directory
    #note: aka is semi-colon separated list of names. check each one.
    cmd = "select id as experiment_id, p.exomekit as kit, p.sample_type as sample_type, p.sample_internal_name as chgvid, p.sample_external_name as origid, p.sample_aka as aka from prepT p, SampleT s, Experiment e where e.id = p.experiment_id and e.sample_id = s.sample_id and p.status like 'External%'"
    samp_info, numrows = run_query_mml_no_commit(cmd, conn_sdb)
    print("there are {} unmapped external chgvid submission".format(numrows))
    
    if numrows == 0:
        print("there are no unnapped external submission ? Check again")
        sys.exit()

    return(samp_info)



def main():

    parser = argparse.ArgumentParser(description=__doc__)    
    parser.add_argument("-m","--manifest",help="optional manifest with optional bam_name column",default="")
    parser.add_argument("-d","--directory_file",help="optioanl filename containing list of directories where to look for files",default="")
    parser.add_argument("-s","--sample_type",help="sample_type",required=True,choices=['exome','genome'])
    args = parser.parse_args() 
    assert os.path.exists(args.directory_file), "file {} exists".format(args.directory_file) 
    setup_logging('/nfs/seqscratch10/mml2204/sample_file_mapping/','sample_file_mapping')
    logger = logging.getLogger(__name__)
    
    dup_chg, dup_orig, dup_aka = [],[],[]
    if args.manifest == "":
        #use DB
        #conn_sdb = get_connection_mml('sequenceDB','sequenceDB')
        #recent_ext_samp = get_recent_ext_sub(conn_sdb)
        
        #conn_sdb.close()
        #code not complete
        pass

    else:
        #parse input manifest 
        #will need to add a "count" in lieu of expt_id
        assert os.path.isfile(args.manifest), "input manifest file {} not found".format(args.manifest)
        
        chg, aka, orig = [], [], []

        with open(args.manifest, newline='') as csvfile:
            manifest_reader = csv.DictReader(csvfile)#delimiter = ',' ? 
            header = manifest_reader.fieldnames
            for entry in imp_fields_from_man:
                assert entry in header, "{} field is not in csv header".format(entry)
            #imp_fields_from_man = ['CHGVID','OrigID','AKA','SeqType','FundCode','SubProject','Protocol','CurrProjLeader','GAFbin','FamilyID','FamilyRelationProband','exomeKit']
            
            """
            first = manifest_reader.next()
            subproject = first['SubProject']
            gaf = first['GAFbin']
            protocol = first['Protocol']
            irb = first['FundCode']
            projleader = first['CurrProjLeader']
            
            chg = [row['CHGVID']]
            aka = [row['AKA']]
            origid = [row['OrigID']]
            """
    

            for row in manifest_reader:
                
                """
                assert row['SubProject'] == subproject, "subproject {} first row subproject: {}".format(row['SubProject'],subproject)
                assert row['GAFbin'] == gaf, "gafbin {} first row gafbin: {}".format(row['GAFbin'],gaf)
                assert row['Protocol'] == protocol, " protocol = {}, first row protocol = {}".format(row['Protocol'],protocol)
                assert row['FundCode'] == irb, "irb = {}, first row irb = {}".format(row['FundCode'],irb)
                assert row['CurrProjLeader'] == projleader, "projleader = {}, first row projleader = {}".format( row['CurrProjLeader'],projleader)
                """

                chg.append(row['CHGVID'])
                orig.append(row['OrigID'])
                aka.append(row['AKA'])
                
                assert row['SeqType'].lower() == args.sample_type.lower()

                if row['SeqType'].lower() == 'exome':
                    assert row['exomeKit'] is not None and row['exomeKit'] != 'N/A', "incorrect exomkit value for {}".format(row['CHGVID'])
                else:
                    assert row['exomeKit'] == 'N/A', "not an exome so should have kit = 'N/A' for {}".format(row['CHGVID'])
                    
                if 'file_name' in row and len(row['file_name']) > 1:
                    print("this has not yet bee coded")
                    pass
                    """
                    filenames = row['file_name'].split('|')
                    for entry in filenames:
                        assert os.path.exists(entry)
                        #get filetype  
                    """
        dup_chg = list(set([x for x in chg if chg.count(x) > 1]))
        #print(dup_chg)
        dup_aka = list(set([x for x in aka if aka.count(x) > 1]))
        #print(dup_aka)
        dup_orig = list(set([x for x in orig if orig.count(x) > 1]))
        #print(dup_orig)

    print("dup chg:{}".format(','.join(dup_chg)))
    print("dup orig:{}".format(','.join(dup_orig)))
    print("dup aka:{}".format(','.join(dup_aka)))
                
    if args.directory_file == "":
        #lookup all dirs in /nfs/tx and then look up redmine to see if the corresponding ticket is open && and has topic "External Transfer Registration"; if not ignore that directory
        print("this has not yet been coded. enter a directory")
        sys.exit()
        #file_dict = get_all_files_in_dir("/nfs/tx")

    else:
        file_dict = get_all_files_in_dir(args.directory_file)
        

    #print("dup chg: {} ".format(','.join(dup_chg)))
    #print("dup orig: {} ".format(','.join(dup_orig)))
    mapping_op = get_file_data_mapping(chg,aka,orig,dup_chg,dup_aka,dup_orig,file_dict) 


main()






