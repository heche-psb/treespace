import pandas as pd
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import logging
from collections import Counter
import os
from joblib import Parallel, delayed
from tqdm import tqdm,trange
import subprocess as sp
from Bio import AlignIO
from Bio.Alphabet import IUPAC

def _mkdir(dirname):
    if not os.path.isdir(dirname) :
        os.mkdir(dirname)
    return dirname

def fetcher(config):
    """
    Return the essential info in dict
    """
    para_dict = {}
    getinfo = lambda x : (x[0],x[1].strip())
    with open(config,'r') as f:
        for line in f.readlines():
            key,value = getinfo(line.split('\t')[0:2])
            para_dict[key] = value
    return para_dict

def reporter(dic):
    """
    Log out the info
    """
    for key,value in dic.items(): logging.info("{0}\t{1}".format(key,value))

def listdir(data,align=False):
    """
    Return the path to sequence files in a list
    """
    parent = os.getcwd()
    os.chdir(data)
    y = lambda x:os.path.join(data,x)
    files_clean = [y(i) for i in os.listdir() if i!="__pycache__"]
    if not align:
        logging.info("In total {} sequence files are found".format(len(files_clean)))
        logging.info(", ".join([os.path.basename(i) for i in files_clean]))
    os.chdir(parent)
    return files_clean

def write_seq(fam,og,seq,seqtype):
    """
    og is a series
    """
    og = og.dropna()
    with open(fam+'.'+seqtype,'w') as f:
        for sp,gids in og.items():
            for gid in gids.split(', '): f.write('>{0}\n{1}\n'.format(gid,seq[sp][gid].seq))

def write_seq_translate(fam,og,seq,to_stop,cds):
    """
    og is a series
    """
    og = og.dropna()
    with open(fam+'.'+'pep','w') as f:
        for sp,gids in og.items():
            for gid in gids.split(', '): f.write('>{0}\n{1}\n'.format(gid,seq[sp][gid].translate(to_stop=to_stop,cds=cds,id=gid).seq))

def find_singleton(og):
    """
    return True or False
    """
    og = og.dropna()
    num = 0
    for seqs in og.values:
        num += len(seqs.split(', '))
        if num > 1:
            return False
    return True

def deal_options(options):
    Options = []
    for i in options:
        if len(i.split(' '))==1: Options.append(i)
        else: Options+=i.split(' ')
    return Options

def mafft(fseq,options,faln):
    if options == '': cmd = ["mafft"] + [fseq]
    else: cmd = ["mafft"] + options + [fseq]
    out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    with open(faln, 'w') as f: f.write(out.stdout.decode('utf-8'))

def muscle(fseq,options,faln):
    if options == '': cmd = ["muscle"] + ['-in',fseq] + ['-out',faln]
    else: cmd = ["muscle"] + ['-in',fseq] + ['-out',faln] + options
    sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)

def prank(fseq,options,faln):
    if options == '': cmd = ["prank"] + ['-d='+fseq] + ['-o='+faln]
    else: cmd = ["prank"] + ['-d='+fseq] + ['-o='+faln] + options
    sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)

def iqtree(faln,options):
    cmd = ["iqtree", "-s", faln] + options
    out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    print(cmd,out.stderr)

def iqtree2(faln,options):
    cmd = ["iqtree2", "-s", faln] + options
    sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)

def fasttree(faln,options):
    cmd = ["FastTree"] + options + ["-out", faln+'.FastTree', faln]
    sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)

def mrbayes(faln,options):
    faln_nex = faln+".nexus"
    AlignIO.convert(faln, 'fasta', faln_nex, 'nexus', IUPAC.extended_protein)
    conf =  faln_nex+".config.mb"
    logf = faln_nex+".mb.log"
    bashf = faln_nex+".bash.mb"
    with open(conf,"w") as f:
        if 'set' in options: f.write('set'+' '+' '.join(options['set'])+'\n')
        else: f.write("set autoclose=yes nowarn=yes\n")
        f.write("execute {}\n".format(os.path.basename(faln_nex)))
        if 'prset' in options: f.write('prset'+' '+' '.join(options['prset'])+'\n')
        else: f.write("prset ratepr=variable\n")
        if 'lset' in options: f.write('lset'+' '+' '.join(options['lset'])+'\n')
        else: f.write("lset rates=gamma\n")
        if 'mcmcp' in options: f.write('mcmcp'+' '+' '.join(options['mcmcp'])+'\n')
        else: f.write("mcmcp diagnfreq=100 samplefreq=10\n")
        if 'mcmc' in options: f.write('mcmc'+' '+' '.join(options['mcmc'])+'\n')
        else: f.write("mcmc ngen=1100 savebrlens=yes nchains=1\n")
        f.write("sumt\nsump\nquit\n")
    with open(bashf,"w") as f:
        f.write('mb <{0}> {1}'.format(conf,logf))
    mb_cmd = ["sh", bashf]
    sp.run(mb_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    with open('rm.sh',"w") as f: f.write("rm *.bash.mb")
    rm_cmd = ['sh', 'rm.sh']
    sp.run(rm_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    rm_cmd = ['rm', 'rm.sh']
    sp.run(rm_cmd, stdout=sp.PIPE, stderr=sp.PIPE)


class Config_Hauler:
    """
    Fetch the info in the config file and implement corresponding analysis
    """
    def __init__(self,config,outdir):
        logging.info("Fetching infomation from config file")
        self.para,self.outdir = fetcher(config),outdir
        self.OG = pd.read_csv(self.para['Orthogroup path:'],header=0,index_col=0,sep='\t')
        self.threads,self.seqtype = int(self.para['Number of threads:']),self.para['Sequences type:']
        self.data = self.para['Sequences directory:']
        reporter(self.para)
        _mkdir(outdir)
        if self.para['Sequences form:'] == 'species':
            self.read_seq()
            logging.info("Writing sequences per family")
            self.write_famseq()
        logging.info("Aligning each family")
        self.aligning()
        logging.info("Inferring gene tree per family")
        self.genetree()

    def read_seq(self):
        seq_paths = listdir(self.data)
        self.gsmap,self.SEQ = {},{}
        for seq in seq_paths:
            self.SEQ[os.path.basename(seq)] = {}
            for record in tqdm(SeqIO.parse(seq, 'fasta'),desc="Reading {}".format(os.path.basename(seq)),unit=" sequences"):
                if not (self.gsmap.get(record.id) is None):
                    logging.error("Duplicated gene id found in {}".format(os.path.basename(seq)))
                    exit(1)
                else:
                    self.gsmap[record.id] = os.path.basename(seq)
                    self.SEQ[os.path.basename(seq)][record.id] = record

    def write_famseq(self):
        OGSEQ_dir = _mkdir(os.path.join(self.outdir,"OG_SEQ"))
        parent = os.getcwd()
        os.chdir(OGSEQ_dir)
        self.Fam_list = Fam_list = list(self.OG.index)
        self.Fam_path = {i:os.path.join(OGSEQ_dir,i+'.'+self.seqtype) for i in Fam_list}
        self.SOGs = list(self.OG.index[[find_singleton(self.OG.loc[fam,:]) for fam in Fam_list]])
        self.NonSOGs = list(set(Fam_list)-set(self.SOGs))
        if len(set(self.OG.index)) != len(Fam_list):
            element_freq = Counter(Fam_list)
            duplicated = [element for element,freq in element_freq.items() if freq != 1]
            logging.error("Duplicated gene family id found for {}".format(", ".join(duplicated)))
            exit(1)
        Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(write_seq)(Fam_list[i],self.OG.loc[Fam_list[i],:],self.SEQ,self.seqtype) for i in trange(len(Fam_list)))
        os.chdir(parent)
        if self.para['Translation:'] == 'yes':
            logging.info("Writing translated sequences per family")
            OGSEQPEP_dir = _mkdir(os.path.join(self.outdir,"OG_SEQ_TRANSLATE"))
            os.chdir(OGSEQPEP_dir)
            to_stop,cds = self.para['Translation parameters:'].split(',') 
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(write_seq_translate)(Fam_list[i],self.OG.loc[Fam_list[i],:],self.SEQ,to_stop,cds) for i in trange(len(Fam_list)))
            self.Fam_trans_path = {i:os.path.join(OGSEQPEP_dir,i+'.pep') for i in Fam_list}
            os.chdir(parent)

    def aligning(self):
        if self.para['Aligner parameters:'] == 'default': options = ''
        else:
            y = lambda x:[i.strip() for i in x]
            options = deal_options(y(self.para['Aligner parameters:'].split(',')))
        self.OGALIGN_dir = OGALIGN_dir = _mkdir(os.path.join(self.outdir,"OG_SEQ_ALIGNMENT_TREE"))
        if self.para['Aligner:'] == 'mafft':
            z = lambda x:os.path.join(OGALIGN_dir,x+".mafft")
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(mafft)(self.Fam_path[self.NonSOGs[i]],options,z(self.NonSOGs[i])) for i in trange(len(self.NonSOGs)))
        elif self.para['Aligner:'] == 'muscle':
            z = lambda x:os.path.join(OGALIGN_dir,x+".muscle")
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(muscle)(self.Fam_path[self.NonSOGs[i]],options,z(self.NonSOGs[i])) for i in trange(len(self.NonSOGs)))
        elif self.para['Aligner:'] == 'prank':
            z = lambda x:os.path.join(OGALIGN_dir,x+".prank")
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(prank)(self.Fam_path[self.NonSOGs[i]],options,z(self.NonSOGs[i])) for i in trange(len(self.NonSOGs)))
        self.Aln_Path = {fam:path for fam,path in zip(sorted(self.NonSOGs),listdir(OGALIGN_dir,align=True))}

    def genetree(self):
        y = lambda x:[i.strip() for i in x]
        if self.para['Tree algorithm parameters:'] == 'default': options = ''
        elif self.para['Tree algorithm:'] == 'mrbayes':
            options = {}
            for p in y(self.para['Tree algorithm parameters:'].split(',')):
                if options.get(p.split(' ')[0]) is None:
                    options[p.split(' ')[0]] = p.split(' ')[1:]
                else:
                    options[p.split(' ')[0]] = options[p.split(' ')[0]]+p.split(' ')[1:]
        else:
            options = deal_options(y(self.para['Tree algorithm parameters:'].split(',')))
        keys = sorted(self.Aln_Path.keys())
        values = {i:self.Aln_Path[j] for i,j in enumerate(keys)}
        if self.para['Tree algorithm:'] == 'iqtree':
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(iqtree)(values[i],options) for i in trange(len(keys)))
        if self.para['Tree algorithm:'] == 'iqtree2':
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(iqtree2)(values[i],options) for i in trange(len(keys)))
        if self.para['Tree algorithm:'] == 'fasttree':
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(fasttree)(values[i],options) for i in trange(len(keys)))
        if self.para['Tree algorithm:'] == 'mrbayes':
            parent = os.getcwd()
            os.chdir(self.OGALIGN_dir)
            Parallel(n_jobs=self.threads,backend='multiprocessing')(delayed(mrbayes)(os.path.basename(values[i]),options) for i in trange(len(keys)))
            os.chdir(parent)





