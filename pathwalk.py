import os
import pdb
import linecache
class pathwalker:
    def __init__(self,pdbname,mapfile,numOfAA,denThr):
        self.pdbname=pdbname
        self.mapfile=mapfile
        self.numOfAA=numOfAA
        self.denThr=denThr
        self.envEMAN="/home/qjq/EMAN2/examples/"
        self.returnValue = returnValue
    def check(self,sheetfile):
        N_sheet = int(linecache.getlines(sheetfile)[-2].split()[-1])+1
        return N_sheet
    def mainAutoProtocal(self):
        start_command = self.envEMAN+"e2pathwalker_auto.py "+ self.mapfile +" --natoms="+ self.numOfAA + " --denthr="+ self.denThr
        os.system(start_command)
            #e2pathwalker_auto.py ../1a6ma.mrc --natoms=150 --denthr=0.471
        start_sheet_dir = "pathwalk_"+"01/"+"sheet.pdb"
        N_input = int(self.numOfAA)
        NAA = int(self.numOfAA)
        N_sheet = check(start_sheet_dir)
        next_command = ""
        justCondition = False
        queryTable=[]
        i = 1
        returnValue = i
        while (not justCondition):
            N_input = N_input + NAA - N_sheet
            next_command = self.envEMAN+ "e2pathwalker_auto.py "+ self.mapfile + " --natoms="+str(N_input)+" --denthr="+ self.denThr
            os.system(next_command)
            if ( i < 10):
                next_sheet_dir = "pathwalk_"+"0"+str(i)+"/sheet.pdb"
            else:
                next_sheet_dir = "pathwalk_"+str(i)+"/sheet.pdb"
            N_sheet = check(next_sheet_dir)
            queryTable.append(abs(N_sheet-NAA))
            if (i <= 20):
                if N_sheet == NAA:
                    justCondition = True
                    returnValue = i
            else:
                returnValue = queryTable.index(min(queryTable))
                justCondition = True
            i = i +1
        return returnValue
