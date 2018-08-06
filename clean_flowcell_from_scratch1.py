import sys
import os
import glob
from utils import *
import warnings



def main():
    if len(sys.argv) != 2:
        raise Exception("need a filename as input argument")

    
    flowcells = []
    with open(sys.argv[1],'r') as f:
        for line in f:
            flowcells.append("'" + line.strip() + "'")
        
    conn_sdb = get_connection_mml('sequenceDB','sequenceDB')    
    cmd = "select * from Flowcell where fcillumid in ({})".format(','.join(flowcells))
    flo_info, num_rows = run_query_mml_no_commit(cmd,conn_sdb)
    assert num_rows  == len(flowcells), "all flowcells in the text file should be in the Flowcell table"

    for entry in flo_info:
        if entry['fc_archive'] is None:
            loc = '/nfs/{}'.format(entry['SeqsataLoc'])
        else:
            loc = '{}'.format(entry['fc_archive'])
      
        cmd = "SELECT DISTINCT CHGVID, lanenum, UPPER(sample_type) as type FROM prepT p JOIN Lane l ON l.PREPID=p.PREPID JOIN Flowcell f ON l.FCID=f.FCID WHERE FCIllumID='{}' and fcillumid not like '%failed%'".format(entry['FCillumID'])
        samp_info, num_rows = run_query_mml_no_commit(cmd,conn_sdb)
        
        for sample in samp_info:
            #fastqs exist
            print(entry['FCillumID'])
            print(sample)
            print(loc)

            if not os.path.exists('{0}/{1}/{2}'.format(loc,sample['type'],sample['CHGVID'])):
                sample['CHGVID'] = sample['CHGVID'] + 'rep'
                warnings.warn("{} This is a rep".format(sample['CHGVID']))

            r1 = glob.glob('{0}/{1}/{2}/{3}/{2}_S*_L00{4}_R1_001.fastq.gz'.format(loc,sample['type'],sample['CHGVID'],entry['FCillumID'],sample['lanenum']))
            r2 = glob.glob('{0}/{1}/{2}/{3}/{2}_S*_L00{4}_R2_001.fastq.gz'.format(loc,sample['type'],sample['CHGVID'],entry['FCillumID'],sample['lanenum']))
            assert len(r1) == 1, "only 1 file required for {}".format(r1[0])
            print(r1)
            assert len(r2) == 1, "only 1 file required for {}".format(r2[0])
            print(r2)
            if not os.path.exists('{0}/{1}/{2}/{3}/laneBarcode.html'.format(loc,sample['type'],sample['CHGVID'],entry['FCillumID'])):

                    warnings.warn("lanebarcode file {0}/{1}/{2}/{3}/laneBarcode.html must exist".format(loc,sample['type'],sample['CHGVID'],entry['FCillumID']))

            #fastq sizes are some minimum
            #sz_r1 = os.path.getsize(r1[0])
            #sz_r2 = os.path.getsize(r2[0])
            
            #assert sz_r1 > 1000000000 and sz_r2 > 1000000000, "{}, {}  must be larger than 1GB".format(r1[0],r2[0])

        
        #check for SAV
        sav_loc = '{}/summary/SAV/'.format(loc)
        sav_file = glob.glob('{}/{}_*_SAV.tar.gz'.format(sav_loc,entry['FCillumID']))
        print(sav_file)
        assert len(sav_file) == 1, "only 1 SAV file per flowcell at {}".format(sav_file[0])


main()
