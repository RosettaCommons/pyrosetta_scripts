import py3Dmol
import pyrosetta.distributed.io
from ipywidgets import interact, IntSlider


__author__ = "Jason C. Klima"
__email__  = "klimaj@uw.edu"


def view_poses(poses=None):
    """View a list of PyRosetta poses in py3Dmol with the ipywidgets interactive slider within a Jupyter notebook.
    The user must have already initialized PyRosetta providing .params files for any ligands in the poses.
    """
    
    def view_py3Dmol(i=0):
        viewer = py3Dmol.view(1200, 800)
        viewer.addModels(pdbstrings[i], "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum"}, "stick": {"radius": 0.2}})
        viewer.setStyle({"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "greenCarbon", "radius": 0.2}})
        viewer.zoomTo()
        return viewer.show()

    if not isinstance(poses, list):
        raise ValueError("Input should be a list of poses.")

    pdbstrings = [pyrosetta.distributed.io.to_pdbstring(p) for p in poses]

    s_widget = IntSlider(min=0, max=len(pdbstrings)-1, description="Structure", continuous_update=False)
    return interact(view_py3Dmol, i=s_widget)


def view_pdbs(pdb_filenames=None):
    """View a list of .pdb files in py3Dmol with the ipywidgets interactive slider within a Jupyter notebook.
    The user must have already initialized PyRosetta providing .params files for any ligands in the .pdb files.
    """

    def view_py3Dmol(i=0):
        viewer = py3Dmol.view(1200, 800)
        viewer.addModels(pdbstrings[i], "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum"}, "stick": {"radius": 0.2}})
        viewer.setStyle({"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "greenCarbon", "radius": 0.2}})
        viewer.zoomTo()
        return viewer.show()

    if not isinstance(pdb_filenames, list):
        raise ValueError("Input should be a list of .pdb file paths.")

    pdbstrings = [pyrosetta.distributed.io.to_pdbstring(pyrosetta.distributed.io.pose_from_file(p)) for p in pdb_filenames]

    s_widget = IntSlider(min=0, max=len(pdbstrings)-1, description="Structure", continuous_update=False)
    return interact(view_py3Dmol, i=s_widget)
