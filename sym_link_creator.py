""" EX: [CHGVID]/raw/cg1096416_FASTQ13.fastq.gz --> [CHGVID]/[fcillumid]/[IGMID]_AAAAAA_L001_R1_001.fastq.gz"""
import argparse
import gzip
import io
import os
import re
import subprocess
import sys
from glob import glob


def main():
    args = parse_arguments()
    sample_path = args.fastq_folder
    sym_loc = args.symlink_folder
    if sym_loc == '':
        sym_loc = sample_path
    verbose = args.verbose

    while  sample_path[-1] == '/':
        sample_path = sample_path[0:-1]
    if os.path.isabs(sample_path) == False:
        raise Exception("Location is not absolute path!")
    sample = sample_path.split('/')[-1]
    if os.path.isdir('{}/raw'.format(sample_path))== False:
        raise Exception("Raw fastq folder not found!")
    check_for_fastqs(sample_path)
    check_for_symlinks(sym_loc,verbose)

    fastqs = glob('{}/raw/*fastq.gz'.format(sample_path))
    if fastqs == []:
        fastqs = glob('{}/raw/*/*fastq.gz'.format(sample_path))
        if fastqs == []:
            fastqs = glob('{}/raw/*txt.gz'.format(sample_path))
            if fastqs == []:
                #hack
                fastqs = glob('{}/raw/*fq*.gz'.format(sample_path))
                if fastqs == []:
                    raise Exception("No fastq.gz were found!")
    total_fastqs = len(fastqs)
    for fastq in fastqs:
        fcillumid,lane,read,adapter = get_rg_info(fastq,verbose)
        fastq_counter = fastq.split('.')[0].split('_')[-1]
        if len(fastq_counter) != 3 or fastq_counter.isdigit() == False:
            fastq_counter = '001'
        sym_link_fastq = ('{}/{}/{}_{}_L00{}_R{}_{}.fastq.gz'
                         ).format(sym_loc,fcillumid,sample,adapter,
                                  lane,read,fastq_counter)
        fastq_dir = '{}/{}'.format(sym_loc,fcillumid)
        make_fcillumid_dir(fastq_dir,verbose)

        ln_cmd = ['ln','-s',fastq,sym_link_fastq]

        if os.path.exists(sym_link_fastq):

            raise Exception('Symlink already exists!: {}'.format(ln_cmd))
        else:
            if verbose:
                print(' '.join(ln_cmd))
            subprocess.check_call(ln_cmd)
    symlinks = glob('{}/*/*fastq.gz'.format(sym_loc))
    total_symlinks = []
    for symlink in symlinks:
        if 'raw' not in symlink:
            total_symlinks.append(symlink)
    total_symlinks_count = len(total_symlinks)
    if total_symlinks_count != total_fastqs:
        print(total_symlinks_count,total_fastqs)
        raise Exception('Number of raw fastqs and symlinks do not match!')
    else:
        print("Sample {}'s symlink creation done".format(sample))

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Display verbose output")
    parser.add_argument('-f','--fastq-folder', required=True,
                        help="Specify scratch dir for bcl2fastq")
    #hack to create symlinks at a different location
    parser.add_argument('-d','--symlink-folder',required=False,default='',
            help='required when symlink is not in same location (base directory) as raw file')
    args=parser.parse_args()
    return args

def check_for_fastqs(sample_path):
    fastqs = glob('{}/raw/*fastq'.format(sample_path))
    if fastqs != []:
        raise Exception("Fastqs not gzip!")

def check_for_symlinks(sample_path,verbose):
    possible_folders = glob('{}/*'.format(sample_path))
    for folder in possible_folders:
        if os.path.isdir(folder) == True and 'raw' not in folder:
            clean_up_symlink_folder(folder,verbose)

def clean_up_symlink_folder(folder,verbose):
    files = glob('{}/*'.format(folder))
    for file in files:
        if os.path.islink(file) == True:
            if verbose:
                print('removing symlink: {}'.format(file))
            os.remove(file)
    files = glob('{}/*'.format(folder))
    if files != []:
        raise Exception("Files still exist in folder")
    else:
        if verbose:
            print('removing folder: {}'.format(folder))
        os.rmdir(folder)


def make_fcillumid_dir(fastq_dir,verbose):
    if not os.path.exists(fastq_dir):
        try:
            if verbose:
                print('mkdir {}'.format(fastq_dir))
            os.makedirs(fastq_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                    raise

def check_fcillumid(fcillumid):
    # As of writing in this script all recent illumina flowcells should begin
    # with C, D, or H and end with X or Y or X2 (for NextSeq High output flowcell).
    # Also the fcillumid should only be 9 characters long.
    search_obj = re.search("^H[A-Z,0-9]{4}B[CG]X2$",fcillumid) #for NextSeq High output
    if search_obj is None:
        search_obj = re.search('[CHD].*[xXyY]',fcillumid)
    if search_obj is None:
        raise Exception('Check fcillumid: {}'.format(fcillumid))
    else:
          fcillumid = search_obj.group()
    if len(fcillumid) != 9:
        raise Exception('Check fcillumid length: {}'.format(fcillumid))
    else:
        return fcillumid

def check_values(machine,fcillumid,lane,tile,adapter,read):
    #Sanity Checks for all values parsed from the read header

    #FCILLUMID check
    check_fcillumid(fcillumid)

    #Tile Check
    # Tiles are structured two ways:
    #
    # NovaSeqs/HiSeqs: [Surface][Swath][Section(1st digit)][Section(2nd digit)]
    # NextSeqs: [surface][swath][camera section][tile digit1][tile digit 2]
    #
    # However determining the machine type via machine name will be very 
    # problematic and unrelible hence I just check the len of the tile
    if tile.isdigit() == False:
        raise Exception('tile is not numeric: {}!'.format(tile))
    if len(tile) == 4: #HiSeq and NovaSeqs
        tile_search_ojb = re.search('^[12][123][012][0-9]$',tile)
    elif len(tile) == 5 and re.match("^@N[BS][0-9]{6}$",machine) is not None: #NextSeqs
        tile_search_ojb = re.match('^[12][123][1-6](0[1-9]|1[0-2])$',tile)
    else:
        raise Exception('Tile is too long!: {}'.format(machine))
    if tile_search_ojb is None:
        raise Exception('Tile is not formatted correctly: {}!'.format(tile))

    #Lane Check
    #Hiseqs should have no more than 2-8 lanes, NovaSeq has 2-4 lanes.
    allowed_lanes =  [str(x) for x in range(1,9)]
    if lane not in allowed_lanes:
        raise Exception('Lane is not formatted correctly: {}!'.format(lane))

    #adapter checks
    allowed_index_bases = ['A','C','T','G','N']
    for base in list(set(adapter)):
        if base not in allowed_index_bases:
            raise Exception('Adapter has illegal bases: {}!'.format(adapter))
    if len(adapter) < 6:
        raise Exception('Adapter is too short: {}!'.format(adapter))

    #read check
    allowed_reads = ['1','2','3']
    if read not in allowed_reads:
        raise Exception('read is not formatted correctly: {}!'.format(read))

def get_rg_info(fastq,verbose):
    gz = gzip.open(fastq,'rt')
    fastq_read_header = gz.readline()
    gz.close()
    rg_info = fastq_read_header.strip().split(':')
    print(rg_info)
    if verbose:
        print(fastq_read_header.strip())
    if len(rg_info) == 10:
        """ IGM Illumina read head example
        @K00347:44:HL7K2BBXX:5:1101:13352:1384 1:N:0:NCTATCCT+NGGATAGG
        @Instrument:RunID:FlowCellID:Lane:Tile:X:Y:UMI ReadNum:FilterFlag:0:IndexSequence or SampleNumber"""
        machine,run_number,fcillumid,lane,tile,X,read,control,something,adapter = rg_info

        read = read.split(' ')[1]
    elif len(fastq_read_header.split(':')) == 7:
        """@HISEQ:549:C6PE0ANXX:2:2305:17233:17109/1
        @Instrument:RunID:FlowCellID:Lane:Tile:X:Y IndexSequence"""
        machine,run_number,fcillumid,lane,tile,X,Y_and_read = rg_info
        read = Y_and_read.split('/')[1]
        adapter = 'AAAAAA'
    elif len(fastq_read_header.split(':')) == 5:
        if re.search('[ACGNT]',rg_info[4]):
            """@FCHNYHLBCXX:1:1207:6982:25608#TGACCAAN/1"""
            fcillumid,lane,tile,X,Y_and_read = rg_info
            machine = None
            fcillumid = fcillumid[3:]
            tmp = Y_and_read.split('#')[1]
            adapter,read = tmp.split('/')
        else:
            fcillumid,lane,tile,X,Y_and_read = rg_info
            fcillumid = fcillumid[3:]
            machine = None
            read = Y_and_read.split('/')[1]
            adapter = 'AAAAAA'
    else:

        print(fastq_read_header)
        raise Exception('Incorrect read header format!')
    adapter = adapter.replace('+','')
    check_values(machine,fcillumid,lane,tile,adapter,read)


    if verbose:
        print(fcillumid,lane,read,adapter)
    #rg_info_sanity_check(fcillumid,lane,read,adapter)
    return fcillumid,lane,read,adapter

if __name__ == '__main__':
    main()
