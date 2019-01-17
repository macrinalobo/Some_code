
# @mml2204
#REMEMBER TO SUPPLY THE REFERENCE FOR CRAM -> BAM

from utils import *
import argparse
import glob

JAVA="/nfs/goldstein/software/jdk1.8.0_05/bin/java"
SAMTOOLS="/nfs/goldstein/software/samtools-1.3/samtools"
PICARD="/nfs/goldstein/software/picard-tools-2.2.1/picard.jar"
PIGZ="/nfs/goldstein/software/pigz-2.3.1/pigz"

"""
POSITION ARGUMENT OPTIONS:
0) Begin
1) Shuffle
2) Revert
3) SamToFastq
4) Pigz
"""

def main(chgvid,temp_folder,input_file,run,ref,samp_type,position):
    
    assert os.path.exists(input_file), "input file {} not found".format(input_file)

    ip_bam_dir = "/".join(input_file.split("/")[0:-1])
    if temp_folder == "": 
        temp_folder = ip_bam_dir + "/BAMCONVERT/" + samp_type.upper()
    else:
        temp_folder = temp_folder + "/" + samp_type.upper()

    print(temp_folder)

    if not os.path.exists(temp_folder):
        cmd = """mkdir -p {}""".format(temp_folder)
        print(cmd)
        if run:
            subprocess_exec(cmd)
    
    ###STEP 0: Cram to bam
    if ref[-3:] == ".gz":
        #check if gzip instead of bgzip, convert to bgzip since samtools accepts only bgzip
        filetype_cmd = """file {}""".format(ref)
        print(filetype_cmd)
        filetype = subprocess_exec(filetype_cmd)
        print(filetype)
        if " gzip " in filetype[0]:
            #convert to bgzip
            bgzip_cmd = """gunzip -c {0} | /nfs/goldstein/software/bin/bgzip > {1}.bgz""".format(ref,ref[0:-3])
            print(bgzip_cmd)
            subprocess_exec(bgzip_cmd)
            ref = ref[0:-3] + ".bgz"

    if input_file[-5:] == ".cram":
        input_bam = "{0}/{1}_cram2bam_op.bam".format(temp_folder,chgvid)
        cram2bam_cmd = """{0} view -bhT {1} {2} -o {3}""".format(SAMTOOLS,ref,input_file,input_bam)
        print(cram2bam_cmd)
        if run and position == 'Begin':
            subprocess_exec(cram2bam_cmd)
            assert os.path.exists(input_bam), "bam file {} not created".format(input_bam)
    else:
        input_bam = input_file
    
    if position == "Begin":
        position="Shuffle"
    
    if position == "Shuffle":
        assert not os.path.exists("{0}/{1}_shuf.bam".format(temp_folder,chgvid)), "shuffle file exists delete it first"
  
    #ref link:https://www.biostars.org/p/8764/ 
    
    ###STEP 1: shuffle
    #/nfs/goldstein/software/samtools-1.3/samtools bamshuf -u
    shuf_cmd = """{0} bamshuf -u {1} {2}/{3}_shuf""".format(SAMTOOLS,input_bam,temp_folder,chgvid)
    print(shuf_cmd)
    if run and position == 'Shuffle':
        subprocess_exec(shuf_cmd)

    #VERIFY LINE COUNTS
    ver_ip_cmd = """{0} view {1} | wc -l""".format(SAMTOOLS,input_bam)
    print(ver_ip_cmd)
    input_bam_wc = 0
    if run:
        input_bam_wc = subprocess_exec(ver_ip_cmd)

    ver_op_cmd = """{0} view {1}/{2}_shuf.bam | wc -l""".format(SAMTOOLS,temp_folder,chgvid)
    print(ver_op_cmd)
    shuf_bam_wc = 0
    if run:
        shuf_bam_wc = subprocess_exec(ver_op_cmd)
        assert input_bam_wc == shuf_bam_wc and input_bam_wc != 0, "shuffle failed;line count before:{}, after:{}".format(input_bam_wc,shuf_bam_wc)    
    
    if position == 'Shuffle':
        position = 'Revert'

    ###STEP 2: RevertSam
    if position == "Revert":
        assert not os.path.exists("{0}/{1}.OQ.sam".format(temp_folder,chgvid)), "OQ file exists delete it first" 
    
    revsam_cmd = """{0} -Xmx16G -jar {1} RevertSam INPUT={2}/{3}_shuf.bam OUTPUT={2}/{3}.OQ.sam TMP_DIR={2} RESTORE_ORIGINAL_QUALITIES=true REMOVE_ALIGNMENT_INFORMATION=true VALIDATION_STRINGENCY=LENIENT >>{2}/{3}_stdout 2>>{2}/{3}_stderr""".format(JAVA,PICARD,temp_folder,chgvid)
    print(revsam_cmd)
    if run and position == 'Revert':
        subprocess_exec(revsam_cmd)

    #VERIFY LINE COUNTS
    ver_op_cmd = """{0} view {1}/{2}.OQ.sam | wc -l""".format(SAMTOOLS,temp_folder,chgvid)
    print(ver_op_cmd)
    rev_sam_wc = 0
    if run:
        rev_sam_wc = subprocess_exec(ver_op_cmd)
        assert shuf_bam_wc == rev_sam_wc, "revert to original qual RevertSam failed; line count before:{}, after:{}".format(shuf_bam_wc,rev_sam_wc)
    
    if position == 'Revert':
        position="SamToFastq"

    ###STEP 3: SamToFastq
    fastqs = glob.glob("{0}/{1}/raw/{1}.fq*".format(temp_folder,chgvid))
    if position == "SamToFastq":
        assert len(fastqs) == 0, "fastqs already exist for this sample: {}. delete them first".format(",".join(fastqs))
    
    create_dir_cmd = "mkdir -p {0}/{1}/raw".format(temp_folder,chgvid)
    print(create_dir_cmd)
    if not os.path.exists("{0}/{1}/raw".format(temp_folder,chgvid)) and run:
        subprocess_exec(create_dir_cmd)

    sam2fq_cmd = """{0} -Xmx128G -jar {1} SamToFastq INPUT={2}/{3}.OQ.sam FASTQ={2}/{3}/raw/{3}.fq1 SECOND_END_FASTQ={2}/{3}/raw/{3}.fq2 INCLUDE_NON_PF_READS=TRUE VALIDATION_STRINGENCY=LENIENT >>{2}/{3}_stdout 2>>{2}/{3}_stderr""".format(JAVA,PICARD,temp_folder,chgvid)
    print(sam2fq_cmd)
    if run and position == "SamToFastq":
        subprocess_exec(sam2fq_cmd)

    fastqs = glob.glob("{0}/{1}/raw/{1}.fq*".format(temp_folder,chgvid))
    if run:
        assert len(fastqs) == 2, "Paired end fastqs incorrectly created:{}".format(','.join(fastqs))

    #VERIFY LINE COUNTS
    fq_lc = 0
    for fastq in fastqs:
        ver_op_cmd = """wc -l {0}""".format(fastq)
        print(ver_op_cmd)
        if run:
            fq1_wc_tmp = subprocess_exec(ver_op_cmd)
            fq1_wc = float(fq1_wc_tmp[0].strip().split()[0])
            assert fq1_wc%4 == 0 and fq1_wc != 0, "incorrect number of lines {} in {}".format(fq1_wc,fastq)
            fq_lc = fq_lc + fq1_wc / 4
    assert int(fq_lc) == int(rev_sam_wc[0]), "SamToFastq failed; line count before:{}, after:{}".format(rev_sam_wc[0],fq_lc)
    
    if position == "SamToFastq":
        position = 'Pigz'

    ###STEP 4: Pigz
    zip_fastqs = glob.glob("{0}/{1}/raw/{1}.fq*.gz".format(temp_folder,chgvid))
    if position == "Pigz":
        assert len(zip_fastqs) == 0, "zipped fastqs already exist for this sample: {}. delete them first".format(",".join(zip_fastqs))
    for fastq in fastqs:
        pigz_cmd = """{0} -v -p 4 {1}""".format(temp_folder,fastq)
        print(pigz_cmd)
        if run and position == "Pigz":
            subprocess_exec(pigz_cmd)
            assert not os.path.exists(fastq) and os.path.exists(fastq+".gz"), "zipped fastq {}.gz incorrectly generated".format(fastq)

    ###STEP 5: Remove intermediate files
    rmshuf_cmd = """rm -v {0}/{1}_shuf.bam""".format(temp_folder,chgvid)
    print(rmshuf_cmd)
    rmOQ_cmd = """rm -v {0}/{1}.OQ.sam""".format(temp_folder,chgvid)
    print(rmOQ_cmd)

    if run:
        subprocess_exec(rmshuf_cmd)
        subprocess_exec(rmOQ_cmd)

    
if __name__ == "__main__":
   
    dir_path = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--chgvid",dest="chgvid",default="",required=True)
    parser.add_argument("-t","--temp_folder",dest='temp_folder',default="",required=False)
    parser.add_argument("-i","--input_bam",dest='input_bam',default="",required=True)
    #parser.add_argument("-s","--symlink_loc",default="/nfs/archive/p2018
    parser.add_argument("-r","--run",dest='run',action='store_true',default=False)
    parser.add_argument("-ref","--ref",dest='ref',default="/scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa",required=False)
    parser.add_argument("-s","--sample_type",dest="samp_type",default="",required=True)
    parser.add_argument("-p","--position",dest="position",default="Begin",required=False,choices=['Begin','Revert','SamToFastq','Shuffle','Pigz'])
    args = parser.parse_args()
    setup_logging(dir_path+'/logs/','bam2fastq_'+args.chgvid)
    main(args.chgvid,args.temp_folder,args.input_bam,args.run,args.ref,args.samp_type,args.position)
    
    logger = logging.getLogger(__name__)
    logger.info("all's well with the world")
    print("all's well with the world")
