from modeller import *

class singleChainModeller :
    def __init__(self, targetSeqPath,targetPdbName):
        self.targetSeqPath = targetSeqPath
        self.targetPdbName = targetPdbName
        self.hdfName = targetPdbName+".hdf5"
        self.buildProfileFile=targetPdbName+'build_profile.prf'
        self.buildProfileFile_ali=targetPdbName+'build_profile.ali'

    def searchTemplate(self):
        log.verbose()
        env = environ()
        sdb = sequence_db(env)
        sdb.convert(seq_database_file = self.targetSeqPath, seq_database_format='PIR',
                    chains_list='ALL', minmax_db_seq_len=(30,4000),
                    clean_sequences=True,outfile= self.hdfName)

    def buildProfile(self):
        log.verbose()
        env = environ()

        sdb = sequence_db(env, seq_database_file=self.hdfName,
                          seq_database_format='BINARY',chain_list='ALL')

        aln = alignment(env)
        aln.append(file= self.targetSequencePath, alignment_format='PIR', align_codes='ALL')

        prf = aln.to_profile()

        prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
                  gap_penalties_1d=(-500, -50), n_prof_iterations=1,
                  check_profile=False, max_aln_evalue=0.01)

        prf.write(file=self.buildProfileFile,profile_format='TEXT')

        aln = prf.to_alignment()
        aln.write(file=self.buildProfileFile_ali,alignment_format='PIR')
