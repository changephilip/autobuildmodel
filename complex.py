import pychimera
import chimera
from chimera import Segger
from Segger import segment_dialog
from Segger import fit_dialog

class functionInChimera:

    def complex_segment(input_pdb_path,input_map_path):
        opened=[]
        opened.append(chimera.openModels.open(input_map_path))
        opened.append(chimera.openModels.open(input_pdb_path))

        d = segment_dialog.()
        d.numSteps.set(0)
        d.Segment()

        x = fit_dialog.Fit_Segments_Dialog()
        x.GroupRegionsByChains()

        CurrentGrouped = x.CurrentSegmentation().grouped_regions()

    def localFit(input_map_path,input_pdb_path):
        chimera.openModels.open(input_map_path)
        chimera.openModels.open(input_pdb_path)
        import FitMap
        from chimera import selection
        from FitMap import fitcmd
        mapid = chimera.openModels.list(id = 0)[0]
        pdbid = chimera.openModels.list(id = 1)[0]

        s1 = selection.ItemizedSelection([pdbid])

        fit_list=FitMap.fitcmd.fitmap(s1,mapid,search=30 ,listFits = False,resolution=20)
        for fit in fit_list:
            va = fit.correlation()
            print va
        return fit_list
    def fitnogui(mapid=1,pdbid=0,search = 30,resolution =20):
        from chimera import openModels as om,selection
        m = om.list(id =mapid)[0]


        
