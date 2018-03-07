import pychimera
import chimera
from chimera import Segger
from Segger import segment_dialog
from Segger import fit_dialog

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
    
