import pymysql
import logging
import os
import subprocess
import sys
from datetime import datetime
import configparser
import traceback
from email.mime.multipart import MIMEMultipart as MM
import smtplib,email,email.encoders,email.mime.text,email.mime.base

def setup_logging(logloc,basename):
    time_stamp='{date:%Y%m%d%H%M%S}'.format( date=datetime.now() )
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s: %(name)s: [%(levelname)s] - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        filename=('{}/{}_{}.log').format(logloc,time_stamp,basename))
    logger = logging.getLogger(__name__)

def send_email(subject,text=''):
    logger = logging.getLogger(__name__)
    emailMsg = MM('alternative')  #email.MIMEMultipart.MIMEMultipart('alternative')
    emailMsg['From'] = 'AAA'
    smtpserver = 'localhost'
    emailMsg['Subject'] = subject
    emailMsg['To'] = ', '.join(['AAA'])
    body='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '
    body +='"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"><html xmlns="http://www.w3.org/1999/xhtml">'
    body +='<body style="font-size:12px;font-family:Verdana"><PRE>'
    body += text
    body += '</PRE></body></html>'
    emailMsg.attach(email.mime.text.MIMEText(body,'html'))
    server = smtplib.SMTP(smtpserver)
    server.sendmail(emailMsg['From'],emailMsg['To'],emailMsg.as_string())
    server.quit()

def get_connection_mml(client,database):
    try:
        reader = configparser.RawConfigParser()
        reader.read('{}/bcl_user.cnf'.format(os.path.dirname(os.path.realpath(__file__))))
        db = 'client' + client
        db_host = reader.get(db, 'host')
        db_user = reader.get(db, 'user')
        db_pass = reader.get(db,'password')
        connection = pymysql.connect(host=db_host,user=db_user,passwd=db_pass,db=database,cursorclass=pymysql.cursors.DictCursor)
        return connection
    except pymysql.Error:
        traceback.print_exc()
        sys.exit("Wrong username/database or password found, please try again- client:{} database:{}".format(client,database))


def run_query_mml(query,connection):
    try:
            cursor = connection.cursor()
            cursor.execute(query)
            results = cursor.fetchall()
            connection.commit()
            return results, cursor.rowcount
    except pymysql.Error as exc:
            raise exc

def run_query_mml_no_commit(query,connection):
    try:
        cursor = connection.cursor()
        cursor.execute(query)
        results = cursor.fetchall()
        return results, cursor.rowcount
    except pymysql.Error as exc:
        raise exc

def subprocess_exec(cmd):
    logger = logging.getLogger(__name__)
    p1 = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    logger.info(cmd)
    out,err = p1.communicate()
    logger.info(out.decode().strip())
    if err.decode() != '':
        logger.info(err.decode())
        raise Exception("problem with {} out:{}, err:{}".format(cmd,out.decode().strip(),err.decode().strip()))
    return out.decode().strip().split('\n')

