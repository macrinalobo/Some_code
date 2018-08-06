import pymysql
import configparser
import logging
import os
import re
import subprocess
import sys
import zipfile
from datetime import datetime
import glob

#checjs that all samples on flowcell are archived to seqsata location. Also rsync with --remvoe-source so fastqs from scratch are removed
#code assumes max of 8 lanes


def setup_logging(logloc):
    time_stamp='{date:%Y%m%d%H%M%S}'.format( date=datetime.now() )
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        filename=('{}/{}_rsync_del.log').format(logloc,time_stamp))
    logger = logging.getLogger(__name__)

def get_connection_mml(database):
    try:
        reader = configparser.RawConfigParser()
        reader.read('{}/bcl_user.cnf'.format(os.path.dirname(os.path.realpath(__file__))))
        db = 'client' + database
        db_host = reader.get(db, 'host')
        db_user = reader.get(db, 'user')
        db_pass = reader.get(db,'password')
        connection = pymysql.connect(host=db_host,user=db_user,passwd=db_pass,db='sequenceDB',cursorclass=pymysql.cursors.DictCursor)
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

def main():
    #optional input file of flowcells to ignore deletion from
    ignore_flowcells = []
    if len(sys.argv) == 2:
        if os.path.exists(sys.argv[1]):
            with open(sys.argv[1],'r') as f:
                for line in f:
                    ignore_flowcells.append(line.strip())
        else:
            raise Exception("file {} doesnt exist".format(sys.argv[1]))
    elif len(sys.argv) == 1 :
        pass
    else:
        raise Exception("too many input args")
    #get loc from config file
    config = configparser.ConfigParser()
    config.read('{}/config.ini'.format(os.path.dirname(os.path.realpath(__file__))))
    loc_config = config.get('locs','bcl2fastq_scratch_drive')
    loc = '/nfs/{}'.format(loc_config)
    logs_dir = config.get('locs','logs_dir')
    flowcells = os.listdir(loc+'/BCL')
    raw_loc = config.get('locs','bcl_dir')
    connection = get_connection_mml('sequenceDB')
    for flo in flowcells:
        flo_name = flo.split('_')[2]
        print(flo_name)
        setup_logging(logs_dir)
        logger = logging.getLogger(flo_name)
        if flo_name in ignore_flowcells:
            logger.info("WARNING!!!!!!! ignoring this flowcell as instructed".format(flo_name))
            continue
        if len(glob.glob('{}/BCL/{}/*/*.fastq.gz'.format(loc,flo))) != 0:
            raise Exception("too many fastqs in {}/BCL/{}".format(loc,flo))
        out = subprocess_exec('find {}/BCL/{}/ -iname "*.gz"'.format(loc,flo))
        for fastq_name in out:
            if 'Undetermined' not in fastq_name:
                raise Exception("{}/BCL/{}/{} should not exist".format(loc,flo,fastq_name))
        #get seqsataloc
        seqsataloc,numrows = run_query_mml("select seqsataloc,machine from Flowcell where fcillumid='{}'".format(flo_name),connection)
        if numrows != 1:
            raise Exception("more than 1 entry for {} in Flowcell table".format(fcillumid))
        check_SAV(flo_name,seqsataloc[0],raw_loc)
        logger.info("SAV check done")
        check_log(flo_name,seqsataloc[0])
        logger.info("nohup file checked")
        check_checkpoints(flo,seqsataloc[0],loc,raw_loc)
        logger.info("passed check_checkpoints")
        check_flowcell_db(flo_name,connection)
        f_rsync,f_tot = check_fastqs(flo,connection,loc)
        logger.info("{} samples archived out of {} for {}. database entries have been fixed".format(f_rsync,f_tot,flo_name))
        if f_rsync == f_tot and f_tot != 0:
            logger.info("cleaning the unaligned bcl dir")
            clean_unaligned_dir(loc,flo,seqsataloc[0])
        else:
            logger.info("Completed work on flowcell {}".format(flo_name))
            logger.info("-----------------------------------------------------------------------")
        #logging.shutdown()
    connection.close()

def clean_unaligned_dir(loc,flo,seqsataloc):
    flo_name = flo.split('_')[2]
    logger = logging.getLogger(__name__)
    cmd="rm -v {0}/BCL/{1}/EmailSent.txt {0}/BCL/{1}/bcl_complete {0}/BCL/{1}/nohup.sge {0}/BCL/{1}/{2}_{3}_BCL.sh".format(loc,flo,flo_name,seqsataloc['seqsataloc'])
    #print(cmd)
    out = subprocess_exec(cmd)
    if os.path.exists('{}/BCL/{}/StorageComplete'.format(loc,flo)):
        cmd="rm -v {}/BCL/{}/StorageComplete".format(loc,flo)
        #print(cmd)
        out = subprocess_exec(cmd)
    cmd="rm -v {}/BCL/{}/Undetermined_S0_L00?_R?_001.fastq.gz".format(loc,flo)
    #print(cmd)
    out = subprocess_exec(cmd)
    rm_undetermined=['{}/BCL/{}/Reports/html/{}/default/Undetermined/all/laneBarcode.html'.format(loc,flo,flo_name),
            '{}/BCL/{}/Reports/html/{}/default/Undetermined/all/lane.html'.format(loc,flo,flo_name),
            '{}/BCL/{}/Reports/html/{}/default/Undetermined/unknown/laneBarcode.html'.format(loc,flo,flo_name),
            '{}/BCL/{}/Reports/html/{}/default/Undetermined/unknown/lane.html'.format(loc,flo,flo_name)]
    cmd="rm -v {}".format(' '.join(rm_undetermined))
    #print(cmd)
    out = subprocess_exec(cmd)
    html_loc='{}/BCL/{}/Reports/html'.format(loc,flo)
    cmd="rm -v {0}/Report.css {0}/index.html {0}/tree.html".format(html_loc)
    #print(cmd)
    out = subprocess_exec(cmd)
    if os.path.exists("/nfs/{}/summary/Stats".format(seqsataloc['seqsataloc'])) is False:
        out = subprocess_exec("mkdir /nfs/{}/summary/Stats".format(seqsataloc['seqsataloc']))
    stat_loc='/nfs/{}/summary/Stats/{}_{}_stats.zip'.format(seqsataloc['seqsataloc'],flo_name,seqsataloc['machine'])
    scratch_stats = '{}/BCL/{}/Stats'.format(loc,flo)
    files_to_zip = ['{}/BCL/{}/LnFractionEmail.txt'.format(loc,flo),
            '{}/DemultiplexingStats.xml'.format(scratch_stats),
            '{}/AdapterTrimming.txt'.format(scratch_stats),
            '{}/ConversionStats.xml'.format(scratch_stats),
            '{}/Stats.json'.format(scratch_stats) ]
    fastq_summary = [ '{}/FastqSummaryF1L{}.txt'.format(scratch_stats,i) for i in range(1,9) if os.path.exists('{}/FastqSummaryF1L{}.txt'.format(scratch_stats,i))] 
    demux_summary = [ '{}/DemuxSummaryF1L{}.txt'.format(scratch_stats,i) for i in range(1,9) if os.path.exists('{}/DemuxSummaryF1L{}.txt'.format(scratch_stats,i))]
    files_to_zip = ' '.join(files_to_zip + fastq_summary + demux_summary)
    cmd="zip -vmT {} {}".format(stat_loc,files_to_zip)
    #print(cmd)
    out = subprocess_exec(cmd)
    if os.path.exists("/nfs/{}/summary/Reports".format(seqsataloc['seqsataloc'])) is False:
        out = subprocess_exec("mkdir /nfs/{}/summary/Reports".format(seqsataloc['seqsataloc']))
    
    files_to_zip_rep = ['{}/BCL/{}/Reports/html/{}/default/all/all/laneBarcode.html'.format(loc,flo,flo_name),
            '{}/BCL/{}/Reports/html/{}/default/all/all/lane.html'.format(loc,flo,flo_name)]
   
    for gaf in next(os.walk('.'))[1]:
        if gaf not in ['Reports','Stats']:
             files_to_zip_rep.append('{}/BCL/{}/Reports/html/{}/{}/all/all'.format(loc,flo,flo_name,gaf))
            #out = subprocess_exec(cmd)
    rep_all = '{}/BCL/{}/Reports/html/{}/all/all/all'.format(loc,flo,flo_name)
    files_to_zip_rep.extend(['{}/laneBarcode.html'.format(rep_all),'{}/lane.html'.format(rep_all)])
    files_to_zip_rep = ' '.join(files_to_zip_rep)
    cmd="zip -vmT /nfs/{0}/summary/Reports/{1}_{2}_all_html.zip {3}".format(seqsataloc['seqsataloc'],flo_name,seqsataloc['machine'],rep_all)
    #print(cmd)
    out = subprocess_exec(cmd)

    cmd="find {}/BCL/{} -type d -empty -print -delete".format(loc,flo)
    #print(cmd)
    out = subprocess_exec(cmd)
    return

def check_flowcell_db(flo_name,connection):
    logger = logging.getLogger(__name__)
    cmd=("SELECT rtaver,hcsver,DateRead1,DateRTA,DateBcl,DateStor,Seqsataloc,CasavaVer,pipelinecomplete,fail,complete,Machine,MachineType,ChemVer "
        "FROM Flowcell where fcillumid='{}'").format(flo_name)
    logger.info(cmd)
    flo_info,numrows = run_query_mml(cmd,connection)
    if numrows != 1:
        raise Exception("error with {}".format(cmd))
    if None in flo_info[0]:
        raise Exception("Flowcell table not updated {}".format(flo_name))
    if '0000-00-00 00:00:00' in flo_info[2:6] :
        raise Exception("Date not updated for {}".format(flo_name))
    if flo_info[0]['pipelinecomplete'] != 1 or flo_info[0]['fail'] != 0 or flo_info[0]['complete'] != 1:
        raise Exception("flowcell {} didnt finish according to database".format(fcillumid))

def check_log(flo_name,seqsataloc):
    logger = logging.getLogger(__name__)
    nohup_loc = '/nfs/{}/summary/bcl_nohup/{}_{}_bcl2fastqLog.zip'.format(seqsataloc['seqsataloc'],flo_name,seqsataloc['machine'])
    logger.info(nohup_loc)
    if os.path.exists(nohup_loc) is False or os.path.getsize(nohup_loc[0]) < 500:
        raise Exception("/nfs/{}/summary/bcl_nohup/{}_{}_bcl2fastqLog.zip doesnt exist".format(seqsataloc[0],flo_name,seqsataloc[1]))
    p1=subprocess.Popen("unzip -c {}".format(nohup_loc),shell=True,stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.PIPE)
    p2=subprocess.Popen("tail -n2",shell=True,stdin=p1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err=p2.communicate()
    out = out.decode().strip()
    logger.info("unzip -c {} | tail -n2, {}".format(nohup_loc,out))
    if 'Processing completed with 0 errors and 0 warnings.' not in out or err.decode() != '':        
        raise Exception("unzip -c {} | tail -n2 gives incorrect output out:{},err:{}".format(nohup_loc,out,err.decode().strip()))
    return

def check_checkpoints(flo,seqsataloc,loc,raw_loc):
    flo_name = flo.split('_')[2]
    logger = logging.getLogger(__name__)
    if not os.path.exists('{}/BCL/{}/bcl_complete'.format(loc,flo)):
        raise Exception("{}/BCL/{}/bcl_complete does not exist".format(loc,flo))
    if not os.path.exists('{}/BCL/{}/EmailSent.txt'.format(loc,flo)):
        raise Exception("{}/BCL/{}/EmailSent.txt does not exist".format(loc,flo))
    if not os.path.exists('{}/BCL/{}/LnFractionEmail.txt'.format(loc,flo)):
        raise Exception("{}/BCL/{}/LnFractionEmail.txt does not exist".format(loc,flo))
    out = subprocess_exec("tail -n1 {}/BCL/{}/LnFractionEmail.txt".format(loc,flo))
    if len(out)!=1 or (out[0] != '<tr><td colspan="7" align="center">&nbsp</td></tr>' and out[0]!='</tr>'):#2nd happens in old format ex: 180203_A00123B_H333TDSXX_Unaligned/LnFractionEmail.txt
        raise Exception('last line of {}/BCL/{}/LnFractionEmail.txt is incorrect. output: {}, expected:<tr><td colspan="7" align="center">&nbsp</td></tr>'.format(loc,flo,out))
    csv_name = glob.glob('/nfs/{}/Sequencing_SampleSheets/{}_*_{}.csv'.format(seqsataloc['seqsataloc'],seqsataloc['machine'],flo_name))
    if len(csv_name) != 1 or  os.path.getsize(csv_name[0]) < 1024: #smalle than 10KB
        raise Exception("Too many csvs/ incorrect size for {} in {}".format(flo_name,seqsataloc['seqsataloc']))
    #check seqscratch1 flowcell exists; in some cases it may have been deleted so this check is useless
    raw_runs = glob.glob('{}/*_{}_*_{}{}'.format(raw_loc,seqsataloc['machine'][0:-1],seqsataloc['machine'][-1],flo_name)) 
    if len(raw_runs) == 0:
        return
    elif len(raw_runs) != 1 or os.path.exists('{}/StorageComplete'.format(raw_runs[0])) is False:
        raise Exception("problem with existence of {}/*_{}_*_{}{}".format(raw_loc,seqsataloc['machine'][0:-1],seqsataloc['machine'][-1],flo_name))
    return

def check_SAV(flo_name,seqsataloc,raw_loc):
    logger = logging.getLogger(__name__)
    print(seqsataloc)
    sav_loc = '/nfs/{}/summary/SAV/{}_{}_SAV.tar.gz'.format(seqsataloc['seqsataloc'],flo_name,seqsataloc['machine'])
    logger.info(sav_loc)
    if os.path.exists(sav_loc) is False or os.path.getsize(sav_loc) < 104857:
        if len(glob.glob('{}/*_{}_*_{}{}'.format(raw_loc,seqsataloc['machine'][0:-1],seqsataloc['machine'][-1],flo_name))) != 0:
            #there is a seqscratch1 to get SAV manually from
            raise Exception("SAV file /nfs/{}/summary/SAV/{}_{}_SAV.tar.gz not proper".format(seqsataloc['seqsataloc'],flo_name,seqsataloc['machine']))
        else:
            logger.info("WARNING!!!!! CANNOT CREATE SAV FOR THIS FLOWCELL SINCE {}/*_{}_*_{}{}) IS DELETED".format(raw_loc,seqsataloc['machine'][0:-1],seqsataloc['machine'][-1],flo_name))
    return

def subprocess_exec(cmd):
    logger = logging.getLogger(__name__)
    p1 = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    logger.info(cmd)
    out,err = p1.communicate()
    logger.info(out.decode().strip())
    if 'rsync' not in cmd and err.decode() != '':
        logger.info(err.decode())
        raise Exception("problem with {} out:{}, err:{}".format(cmd,out.decode().strip(),err.decode().strip()))
    if 'rsync' in cmd and err.decode() != '':
        logger.info(err.decode())
    return out.decode().strip().split('\n')

#it seems sometimes the same fastq name is generated twice by the DB due to DB data problem. Hence the need for set of all fastq names.
def check_fastqs(flo,connection,loc):
    flo_name=flo.split('_')[2]
    logger = logging.getLogger(__name__)
    cmd=("SELECT distinct upper(p.sample_type) sample_type, p.chgvid,lanenum,seqsataloc,status,rg_status,step_status,p.prepid,failR1,failR2,p_prepid,l.fcid,l.seqid,l.poolid,s.GAFbin gaf,p.adapterLet adapter,p.is_released "
         "FROM prepT p, Lane l, Flowcell f, SampleT s, Experiment e "
         "WHERE f.fcid=l.fcid and l.prepid=p.prepid and f.fcillumid='{}' and s.sample_id=e.sample_id and e.id=p.experiment_id").format(flo_name) 
    logger.info(cmd)
    fastq_names = set()
    chgs,num_chg = run_query_mml(cmd,connection)
    logger.info("there are {} samples on flowcell {}".format(num_chg,flo_name))
    if num_chg == 0:
        raise Exception("Flowcell {} has 0 samples".format(flo_name))
    chgvids = []
    num_fastqs_rsynced = 0
    num_fastqs_a = {}
    num_fastqs_s = {}
    num_fastqs_db = {}
    tot_fastqs_a = 0

    for item in chgs:
        #print(item)
        chgvids.append((item['sample_type'],item['chgvid'],item['lanenum']))
        f1 = glob.glob('{0}/{1}/{2}/{3}/{2}_S*_L00{4}_R1_001.fastq.gz'.format(loc,item['sample_type'],item['chgvid'],flo_name,item['lanenum']))
        f2 = glob.glob('{0}/{1}/{2}/{3}/{2}_S*_L00{4}_R2_001.fastq.gz'.format(loc,item['sample_type'],item['chgvid'],flo_name,item['lanenum']))
        #print(f1[0],f2[0])
        if len(f1) != 1 or len(f2) != 1:
            raise Exception("incorrect number of fastqs in {0}/{1}/{2}/{3}/{2}_".format(loc,item['sample_type'],item['chgvid'],flo_name))
        f1_a = glob.glob('/nfs/{0}/{1}/{2}/{3}/{2}_S*_L00{4}_R1_001.fastq.gz'.format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name,item['lanenum']))
        f2_a = glob.glob('/nfs/{0}/{1}/{2}/{3}/{2}_S*_L00{4}_R2_001.fastq.gz'.format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name,item['lanenum']))
        #print(f1_a,f1_b)
        if len(f1_a) != 1 or len(f2_a) != 1:
            raise Exception("incorrect number of fastqs in /nfs/{0}/{1}/{2}/{3}/{2}_".format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name))
        #checksize
        if os.path.getsize(f1[0]) == 0 or  os.path.getsize(f1[0]) !=  os.path.getsize(f1_a[0]):
            raise Exception("{} size:{}, {} size:{}".format(f1[0],os.path.getsize(f1[0]),f1_a[0],os.path.getsize(f1_a[0])))
        #print(f1[0],f1_a[0],f2[0],f2_a[0])
        if item['failR1'] is not None or item['failR2'] is not None:
            raise Exception("Flowcell {} failed in lane {}".format(flo_name,item['lanenum'])) 
        
        #add this in for final run
        out = subprocess_exec("gzip -t {} {}".format(f1[0],f2[0]))

        if item['status'] in ['Archiving','BCL Started','BCL Complete','BCL']:
            cmd="UPDATE prepT set status='Storage',status_time=unix_timestamp() where prepid='{}'".format(item['prepid'])
            logger.info("current status is {}. cmd is {}".format(item['status'],cmd))
            
            _,aff_rows = run_query_mml(cmd,connection)
            if aff_rows != 1:
                raise Exception("{} returned {} rows".format(cmd,aff_rows))
            
            cmd2=("INSERT INTO statusT (STATUS_TIME,STATUS,PREPID,USERID,POOLID,SEQID) "
                  "SELECT DISTINCT UNIX_TIMESTAMP(),'Storage',pt.PREPID,'16940058',0,0 "
                  "FROM prepT pt WHERE prepid='{}'").format(item['prepid'])
            logger.info(cmd2)
            
            _,aff_rows = run_query_mml(cmd2,connection)
            if aff_rows != 1:
                raise Exception("{} returned {} rows".format(cmd2,aff_rows))
            
        if item['rg_status'] in ['sequencing','fastq_ready'] or item['step_status'] != 'in storage':
            cmd=("UPDATE Lane set rg_status='fastq_copied', step_status='in storage' where "
                "prepid='{}' and lanenum='{}' and fcid='{}' and seqid='{}' and  poolid='{}'").format(item['prepid'],item['lanenum'],item['fcid'],item['seqid'],item['poolid'])
            logger.info("current status is ({},{}). cmd is {}".format(item['rg_status'],item['step_status'],cmd))
            
            _,aff_rows = run_query_mml(cmd,connection)
            if aff_rows != 1:
                raise Exception("{} updated too many {} rows".format(cmd,aff_rows))
            
        
        if item['chgvid'] not in num_fastqs_s:
            num_fastqs_s[item['chgvid']] = len(glob.glob('{}/{}/{}/{}/*.fastq.gz'.format(loc,item['sample_type'],item['chgvid'],flo_name)))
            num_fastqs_a[item['chgvid']] = (len(glob.glob('/nfs/{}/{}/{}/{}/*.fastq.gz'.format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name))),item['sample_type'])
            num_fastqs_db[item['chgvid']] = 2
            tot_fastqs_a += num_fastqs_s[item['chgvid']]
            if num_fastqs_a[item['chgvid']][0] != num_fastqs_s[item['chgvid']]:
                raise Exception("incorrect # fastqs in {0}/{1}/{2}/{3}, /nfs/{4}/{1}/{2}/{3} ({5},{6})".format(loc,item['sample_type'],item['chgvid'],flo_name,item['seqsataloc'],num_fastqs_s[item[1]],num_fastqs_a[item[1]][0]))
            gaf = item['gaf'].replace(" ","") #gaf bin can have space between it's characters but dirname doesnt
            path_scratch = '{}/BCL/{}/Reports/html/{}/{}/{}/{}'.format(loc,flo,flo_name,gaf,item['chgvid'],item['adapter'])
            path_dest = '/nfs/{}/{}/{}/{}/lane_{}.zip'.format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name,item['adapter'])
            
            if os.path.exists('{}/laneBarcode.html'.format(path_scratch)) and os.path.exists('{}/lane.html'.format(path_scratch)):
                cmd="zip -vmT {0} {1}/laneBarcode.html {1}/lane.html".format(path_dest,path_scratch)
                #print(cmd)
                out = subprocess_exec(cmd)
            else:
                raise Exception("laneBarcode / lane file not found in {}".format(path_scratch))

            path_scratch = '{}/BCL/{}/Reports/html/{}/{}/{}/all/laneBarcode.html'.format(loc,flo,flo_name,gaf,item['chgvid'])
            path_ar_2 = '/nfs/{}/{}/{}/{}/laneBarcode.html'.format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name)
            path_ar_2_lane = '/nfs/{}/{}/{}/{}/lane.html'.format(item['seqsataloc'],item['sample_type'],item['chgvid'],flo_name)
            path_ar_1 = '{}/{}/{}/{}/laneBarcode.html'.format(loc,item['sample_type'],item['chgvid'],flo_name)
            
            if (os.path.exists(path_scratch) and os.path.exists(path_ar_2) and os.path.exists(path_ar_1) and os.path.getsize(path_ar_2) > 100 
                and os.path.getsize(path_scratch) == os.path.getsize(path_ar_2) and os.path.getsize(path_ar_2) == os.path.getsize(path_ar_1)):
                    path_scratch_lane = '{}/BCL/{}/Reports/html/{}/{}/{}/all/lane.html'.format(loc,flo,flo_name,gaf,item['chgvid'])
                    cmd="rsync -azpvv --remove-source-files -L --no-owner --no-perms {} {}".format(path_scratch,path_ar_2)
                    out = subprocess_exec(cmd)
                    cmd="rsync -azpvv --remove-source-files -L --no-owner --no-perms {} {}".format(path_ar_1,path_ar_2)
                    out = subprocess_exec(cmd)
                    if os.path.exists(path_scratch_lane):
                        cmd="rsync -azpvv --remove-source-files -L --no-owner --no-perms {} {}".format(path_scratch_lane,path_ar_2_lane)
                        out = subprocess_exec(cmd)
                    else:
                        raise Exception("{} does not exist".format(path_scratch_lane))
            else:
                    raise Exception("problem with laneBarcode.html files {} {} {}".format(path_scratch,path_ar_2,path_ar_1))
        
        else:
            num_fastqs_db[item['chgvid']] += 2
        
        fastq_names.update([f1[0],f2[0]])
        cmd="select is_merged from dragen_sample_metadata where pseudo_prepid='{}'".format(item['p_prepid'])
        logger.info(cmd)
        is_merged,aff_rows = run_query_mml(cmd,connection)
        if aff_rows > 1:
            raise Exception("cmd {} returned {} rows".format(cmd,aff_rows))
        
        if item['sample_type'] == 'RNASEQ' and item['is_released'] < 3:
            logger.info("Not completed dragen alignment")
            logger.info("{0}\t{1}\t{2}\t/nfs/{3}/{0}/{1}/{4}\t{5}/{0}/{1}/{4}\tis_released:{6}\n".format(item['sample_type'],item['chgvid'],item['status'],item['seqsataloc'],flo_name,loc,item['is_released']))
            continue
        elif len(is_merged) == 0 or is_merged[0]['is_merged'] <= 1 or  (is_merged[0]['is_merged'] >= 80000 and is_merged[0]['is_merged'] <= 80100):
            logger.info("Not completed dragen alignment")
            logger.info("{0}\t{1}\t{2}\t/nfs/{3}/{0}/{1}/{4}\t{5}/{0}/{1}/{4}\n".format(item['sample_type'],item['chgvid'],item['status'],item['seqsataloc'],flo_name,loc))
            continue
        else:
            pass
        
        cmd="rsync -azpvv --remove-source-files -L --no-owner --no-perms {} {}".format(f1[0],f1_a[0])
        cmd2="rsync -azpvv --remove-source-files -L --no-owner --no-perms {} {}".format(f2[0],f2_a[0])
        out = subprocess_exec(cmd)
        out = subprocess_exec(cmd2)
        num_fastqs_rsynced += 2
        
    if len(fastq_names) != tot_fastqs_a:
        raise Exception("fastqs names built from DB not same as fastqs present in archive/scratch")
        
    for item in num_fastqs_s:
        if  num_fastqs_a[item][0] != num_fastqs_db[item]:
            logger.info("WARNING!!! incorrect number of fastqs archive:{} according to DB:{} for {}. Problem with DB data or with fastqs in folders?".format(num_fastqs_a[item][0],num_fastqs_db[item],item))
        cmd="rmdir -v --ignore-fail-on-non-empty {}/{}/{}/{}".format(loc,num_fastqs_a[item][1],item,flo_name)
        #print(cmd)
        out = subprocess_exec(cmd)

    logger.info("{} fastqs were rsynced out of {}".format(num_fastqs_rsynced,len(fastq_names))) 
    logger.info('\n'.join(fastq_names))    
    return (num_fastqs_rsynced,len(fastq_names))

main()
