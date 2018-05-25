from modeller import *
import os
class singleChainModeller :
    def __init__(self, targetSeqPath,targetPdbName,numTemplate=7,mapfile,resolution,voxelSize):
        self.targetSeqPath = targetSeqPath
        self.targetPdbName = targetPdbName
        self.numTemplate = numTemplate
        self.ModellerDatabase = "pdb_95.pir"
        self.hdfName = targetPdbName+".hdf5"
        self.buildProfileFile=targetPdbName+'build_profile.prf'
        self.buildProfileFile_ali=targetPdbName+'build_profile.ali'
        self.selectPdbName = ""
        self.selectPdbChain = ""
        self.alignFile_ali = ""
        self.alignFile_pap = ""
        self.mapfile = mapfile
        self.resolution = resolution
        self.voxelSize = voxleSize#A per pixel
    def searchTemplate(self):
        log.verbose()
        env = environ()
        sdb = sequence_db(env)
        if not os.path.exists(self.ModellerDatabase):
            print("ERROR! No pdb_95.pir\n")
            print("Goto salilab.org/modeller/supplemental.html and download it please\n")
        sdb.convert(seq_database_file =self.ModellerDatabase , seq_database_format='PIR',
                    chains_list='ALL', minmax_db_seq_len=(30,4000),
                    clean_sequences=True,outfile= self.hdfName)

    def buildProfile(self):
        log.verbose()
        env = environ()

        sdb = sequence_db(env, seq_database_file=self.hdfName,
                          seq_database_format='BINARY',chains_list='ALL')

        aln = alignment(env)
        aln.append(file= self.targetSeqPath, alignment_format='PIR', align_codes='ALL')

        prf = aln.to_profile()

        prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
                  gap_penalties_1d=(-500, -50), n_prof_iterations=1,
                  check_profile=False, max_aln_evalue=0.01)

        prf.write(file=self.buildProfileFile,profile_format='TEXT')

        aln = prf.to_alignment()
        aln.write(file=self.buildProfileFile_ali,alignment_format='PIR')

    def parsePrf(self):
        #import linecache
        if os.path.exists(self.buildProfileFile):
            FileToParse = open(self.buildProfileFile)
        pdbNameList=[]
        chainNumList=[]
        scoreList=[]
        prfList=[]
        for line in FileToParse:
            if line!=[]:
                if line[0]!="#":
                    linet=line.split()
                    pdbName=linet[1][0:4]
                    chainNum=linet[1][-1]
                    score=int(linet[10].split('.')[0])
                    prfList.append([pdbName,chainNum,score])
        prfList.sort(key= lambda arg : arg[2] )
        TopList=prfList[-self.numTemplate:]
        FileToParse.close()
        return TopList

    def __getSinglePdb(self,PDBname):
        headers = {
            "User-Agent": "pMozilla/5.0 (Macintosh; Intel Mac OS X 10_12_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/61.0.3163.100 Safari/537.36"
}
        import socket
        import urllib
        import requests
        socket.setdefaulttimeout(180)
        url_pre = 'http://files.rcsb.org/download/'
        pdb_name = PDBname + '.pdb'
        cif_name = PDBname + '.cif'
        print("Download "+PDBname+'\n')
        reponse = requests.get(url_pre+pdb_name,headers = headers,timeout = 10).status_code
        if reponse == 200:
             try:
                 urllib.urlretrieve(url_pre+pdb_name,pdb_name)
             except urllib.ContentTooShortError,socket.timeout:
                 print 'Network conditions is not good.Reloading.'
                 print PDBname
                 getpdbsingle(PDBname)
        else :
             if requests.get(url_pre+cif_name,headers = headers,timeout = 10).status_code ==200:
                 try:
                     urllib.urlretrieve(url_pre+cif_name,cif_name)
                 except urllib.ContentTooShortError,socket.timeout:
                     print 'Network conditions is not good.Reloading.'
                     print PDBname
                     getSinglePdb(PDBname)
        

    def compare(self):
        #Download atom files
        if os.path.exists("atom_files")==False:
            os.mkdir("atom_files")
        downloadList = self.parsePrf()
        print downloadList

        os.chdir("atom_files")
        for item in downloadList:
            self.__getSinglePdb(item[0])
        os.chdir("..")

        env = environ()
        env.io.atom_files_directory = [".","atom_files"]

        aln = alignment(env)
        for (pdb,chain,score) in downloadList:
            m = model(env, file=pdb, model_segment=('FIRST:'+chain,'LAST:'+chain))
            aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)

        aln.malign()
        aln.malign3d()
        aln.compare_structures()
        aln.id_table(matrix_file='family.mat')
        env.dendrogram(matrix_file='family.mat',cluster_cut=-1.0)

    def selectTemplate(self):
        self.selectPdbName=""
        self.selectPdbChain=""
        self.alignFile_ali = self.targetPdbName+"-"+self.selectPdbName+self.selectPdbChain+".ali"
        self.alignFile_ali = self.targetPdbName+"-"+self.selectPdbName+self.selectPdbChain+".pap"
        return

    def alignWithStructure(self):
        env = environ()
        env.io.atom_files_directory = ['.','atom_files']

        aln = alignment(env)

        mdl = model(env, file = self.selectedPdb[0], model_segment = ('FIRST:'+self.selectPDB[1],'LAST:'+self.selectPdb[1]))
        aln.append_model(md1, align_codes = self.selectPdb[0]+self.selectPdb[1], atom_files= self.selectPdb[0])

        aln.append(file=self.targetSeqPath,align_codes=self.targetPdbName)

        aln.align2d(max_gap_length=40)

        aln.write(file=self.alignFile_ali,alignment_format='PIR')
        aln.write(file=self.alignFile_pap,alignment_format='PAP')

    def buildModel(self):
        from modeller.automodel import *

        env = environ()
        env.io.atom_files_directory=['.',"atom_files"]
        a = automodel(env, alnfile = self.alignFile_ali,
                      knowns=self.selectPdbName+self.selectPdbChain,
                      sequence=self.targetPdbName,
                      assess_methods=assess.DOPE)
        a.starting_mode = 1
        a.ending_model =5
        a.make()

    def simpleFit(self):
        import mrcfile
        log.verbose()
        env = environ()

        #get arguments
        boxSize = 48


import unittest

class Testunit(unittest.TestCase):

    def test_init(self):

        return
