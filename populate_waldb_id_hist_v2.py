from utils import *
import datetime
import json
from email.mime.multipart import MIMEMultipart as MM
import smtplib,email,email.encoders,email.mime.text,email.mime.base

#the optional input arg is a text file containing 1 pseudo_prepid per line
def main():
    f_supersede = []
    if len(sys.argv) == 1:
        pass
    elif len(sys.argv) == 2:
        with open(sys.argv[1],'r') as f:
            for line in f:
                f_supersede.append(line.strip())
    else:
        raise Exception("incorrect number of input args {}".format(len(sys.argv)))
    conn_waldb = get_connection_mml('waldbm','WalDB')
    conn_sdb = get_connection_mml('sequenceDB','sequenceDB')
    setup_logging('/nfs/seqscratch10/mml2204/DEPRECATION/','waldb_hist_field_update')
    logger = logging.getLogger(__name__)
    email_file = []
    email_file = populate_waldb_id_hist(conn_waldb,conn_sdb,f_supersede,email_file)
    if len(email_file) == 0:
        send_email("No harmful discrepancies between waldb finished = 1 and failure = 0 samples and dsm")
    else:
        send_email("PROBLEM WITH WALDB DATA!",''.join(list(set(email_file))))
    conn_waldb.close()
    conn_sdb.close()

def populate_waldb_id_hist(conn_waldb,conn_sdb,f_supersede,email_file):    
    logger = logging.getLogger(__name__)
    cmd="SELECT sample_id,sample_name,UPPER(sample_type) as sample_type,UPPER(capture_kit) as capture_kit,prep_id,sample_finished,sample_failure from sample"
    logger.info(cmd)
    cursor = conn_waldb.cursor()
    cursor.execute(cmd)
    row=cursor.fetchone()
    if row is None:
        raise Exception("no rows in WalDB sample table ??")
    count = 0
    while row is not None:
        cmd="SELECT sample_name,UPPER(sample_type) as sample_type,UPPER(capture_kit) as capture_kit,pseudo_prepid,waldb_id_hist,is_merged from dragen_sample_metadata where pseudo_prepid={}".format(abs(row['prep_id']))
        logger.info(cmd)
        dsm_op,rows = run_query_mml(cmd,conn_sdb)
        logger.info(dsm_op)
        logger.info(row)
        if row['sample_finished'] == 1 and row['sample_failure'] != 0 and  rows != 0 and dsm_op[0]['is_merged'] in [40,100]:
            print('sample_name={};pp={};waldb_sample_type={};waldb_kit={}'.format(row['sample_name'],row['prep_id'],row['sample_type'],row['capture_kit'] ))
            msg = 'Seems complete in dsm  (is_merged = {} ) but sample_failure:{}, sample_finished:{} in waldb for sample_name={};prep_id={};waldb_sample_type={}waldb_kit={}\n'. \
                    format(dsm_op[0]['is_merged'],row['sample_failure'],row['sample_finished'],row['sample_name'],row['prep_id'],row['sample_type'],row['capture_kit'] )
            email_file.append(msg)
        if rows == 0:
            if row['prep_id'] < 0:
                logger.info('sample was deprecated, deleted from dsm and not rerun {} {} {} {}'.format(row['sample_name'],row['sample_type'],row['capture_kit'],abs(row['prep_id'])))
                #email_file.append('sample was deprecated, deleted from dsm and not rerun {} {} {} {}\n'.format(row['sample_name'],row['sample_type'],row['capture_kit'],abs(row['prep_id'])))
                row=cursor.fetchone()
                continue
            elif row['sample_finished'] == 1 and row['sample_failure'] == 0:
                print("not found in dsm table:sample_name={};pp={};waldb_sample_type={};waldb_kit={}".format(row['sample_name'],row['prep_id'], row['sample_type'],row['capture_kit']))
                email_file.append("SERIOUS!!!! since finished = 1 and failure = 0 but absent in dsm not found in dsm table:sample_name={};pp={};waldb_sample_type={};waldb_kit={}\n".format(row['sample_name'],row['prep_id'], row['sample_type'],row['capture_kit']))
                row=cursor.fetchone()
                continue
            else:
                logger.info("not found in dsm table:sample_name={};pp={};waldb_sample_type={};waldb_kit={}".format(row['sample_name'],row['prep_id'], row['sample_type'],row['capture_kit']))
                print("not found in dsm table:sample_name={};pp={};waldb_sample_type={};waldb_kit={}".format(row['sample_name'],row['prep_id'], row['sample_type'],row['capture_kit']))
                #email_file.append("not found in dsm table:sample_name={};pp={};waldb_sample_type={};waldb_kit={}\n".format(row['sample_name'],row['prep_id'], row['sample_type'],row['capture_kit']))
                row=cursor.fetchone()
                continue
        elif rows > 1:
            raise Exception("Too many rows:{} returned for {}".format(rows,cmd))
        dsm_op = dsm_op[0]
        if dsm_op['sample_name'] != row['sample_name'] or dsm_op['sample_type'] != row['sample_type'] or dsm_op['capture_kit'] != row['capture_kit']:
            #deal with case when old genome now released as GENOME_AS_FAKE_EXOME
            #print(row['prep_id'],dsm_op['sample_name'],row['sample_name'],row['sample_type'],dsm_op['sample_type'],row['capture_kit'],row['capture_kit'],dsm_op['capture_kit'])
            if row['prep_id'] < 0 and dsm_op['sample_name'] == row['sample_name'] and row['sample_type'] == 'GENOME' and \
              dsm_op['sample_type'] == 'GENOME_AS_FAKE_EXOME' and (row['capture_kit'] is None or row['capture_kit'] == 'N/A' or row['capture_kit']=='') and dsm_op['capture_kit'] == 'ROCHE':
                pass
            elif dsm_op['sample_name'] == row['sample_name'] and ((row['sample_type'] == 'GENOME_AS_FAKE_EXOME' and \
              dsm_op['sample_type'] == 'GENOME') or (row['sample_type'] == 'GENOME' and dsm_op['sample_type'] == 'GENOME_AS_FAKE_EXOME')):
                print('chgvid={};pp={};dsm_sample_type={};waldb_sample_type={};dsm_kit={};waldb_kit={}'. \
                        format(dsm_op['sample_name'],row['prep_id'],dsm_op['sample_type'],row['sample_type'],dsm_op['capture_kit'],row['capture_kit'] ))
            else:
                email_file.append("{} {} {} {} {} {} data doesnt match in dsm, waldb sample\n" \
                        .format(dsm_op['sample_name'], row['sample_name'],dsm_op['sample_type'],row['sample_type'],dsm_op['capture_kit'],row['capture_kit']))
                raise Exception("{} {} {} {} {} {} data doesnt match in dsm, waldb sample"
                    .format(dsm_op['sample_name'], row['sample_name'],dsm_op['sample_type'],row['sample_type'],dsm_op['capture_kit'],row['capture_kit']))
        
        #create entry for waldb_id_hist
        if row['prep_id'] > 0 and row['sample_finished'] == 0 and row['sample_failure'] == 0:
            cur_state = {'sample_id':row['sample_id'],'date':str(datetime.datetime.now()),'status':'Weird'}
        elif row['prep_id'] < 0 or row['sample_finished'] != 1 or row['sample_failure'] != 0:
            cur_state = {'sample_id':row['sample_id'],'date':str(datetime.datetime.now()),'status':"Deprecated"}
        elif row['prep_id'] in f_supersede:
            cur_state= {'sample_id':row['sample_id'],'date':str(datetime.datetime.now()),'status':"Superseded"}
        else:
            cur_state= {'sample_id':row['sample_id'],'date':str(datetime.datetime.now()),'status':"Current"}
            
        logger.info('{} {}'.format(dsm_op['waldb_id_hist'],cur_state))
        #check last entry
        
        if dsm_op['waldb_id_hist'] is None:
            #raise Exception("shouldnt be here")
            cmd="UPDATE dragen_sample_metadata set waldb_id_hist='{}' WHERE pseudo_prepid={}".format(json.dumps([cur_state]),dsm_op['pseudo_prepid'])
            #print(cmd)
            logger.info(cmd)
            _,aff_rows = run_query_mml(cmd,conn_sdb)
            if aff_rows != 1:
                raise Exception("Problem executing query {}".format(cmd))
        else:
            curr_db_state = json.loads(dsm_op['waldb_id_hist'])
            count_state, state_exists, seen_sid = 0,0,0
            logger.info(cur_state)
            
            while count_state < len(curr_db_state):
                #print(row['prep_id'])
                if curr_db_state[count_state]['sample_id'] == cur_state['sample_id'] and curr_db_state[count_state]['status'] == cur_state['status']:
                    logger.info("waldb_id_hist field is already updated for pseudo_prepid {}, sample_id {}".format(dsm_op['pseudo_prepid'],cur_state['sample_id']))
                    state_exists += 1
                    seen_sid += 1
                elif curr_db_state[count_state]['sample_id'] == cur_state['sample_id'] and curr_db_state[count_state]['status'] != cur_state['status']:
                    seen_sid += 1
                    curr_db_state[count_state]['status'] = cur_state['status']
                    curr_db_state[count_state]['date'] = str(datetime.datetime.now())

                count_state += 1

            if state_exists > 1 or seen_sid > 1:
                raise Exception("sample_id occurs too many times in waldb_id_hist for pp {}".format(dsm_op['pseudo_prepid']))

            if seen_sid == 0 and state_exists == 0:
                curr_db_state.append(cur_state)
            if state_exists == 0:
                cmd="UPDATE dragen_sample_metadata set waldb_id_hist='{}' WHERE pseudo_prepid={}".format(json.dumps(curr_db_state),dsm_op['pseudo_prepid'])
                logger.info(cmd)
                #print(cmd)
                _,aff_rows = run_query_mml(cmd,conn_sdb)
                if aff_rows != 1:
                    raise Exception("Problem executing query {}".format(cmd))
        row=cursor.fetchone()
        logger.info(dsm_op['pseudo_prepid'])
    return email_file

main()
