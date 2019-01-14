import py3Dmol
import pyrosetta
import pyrosetta.distributed.io
from ipywidgets import interact, IntSlider


__author__ = "Jason C. Klima"
__email__  = "klimaj@uw.edu"


def view_poses(poses=None):
    """View a list of PyRosetta poses in py3Dmol with the ipywidgets interactive slider within a Jupyter notebook.
    The user must have already initialized PyRosetta providing .params files for any ligands/non-canonical residues
    in the poses.
    """
    
    def view_py3Dmol(i=0):
        viewer = py3Dmol.view(1200, 800)
        viewer.addModels(pdbstrings[i], "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum"}, "stick": {"radius": 0.25}})
        viewer.setStyle({"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}})
        viewer.zoomTo()
        return viewer.show()

    if not isinstance(poses, list):
        raise ValueError("Input should be a list of poses.")

    pdbstrings = [pyrosetta.distributed.io.to_pdbstring(p) for p in poses]

    s_widget = IntSlider(min=0, max=len(pdbstrings)-1, description="Structure", continuous_update=False)
    return interact(view_py3Dmol, i=s_widget)


def view_pdbs(pdb_filenames=None):
    """View a list of .pdb files in py3Dmol with the ipywidgets interactive slider within a Jupyter notebook.
    The user must have already initialized PyRosetta providing .params files for any ligands/non-canonical residues
    in the .pdb files.
    """

    def view_py3Dmol(i=0):
        viewer = py3Dmol.view(1200, 800)
        viewer.addModels(pdbstrings[i], "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum"}, "stick": {"radius": 0.25}})
        viewer.setStyle({"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}})
        viewer.zoomTo()
        return viewer.show()

    if not isinstance(pdb_filenames, list):
        raise ValueError("Input should be a list of .pdb file paths.")

    pdbstrings = [pyrosetta.distributed.io.to_pdbstring(pyrosetta.distributed.io.pose_from_file(p)) for p in pdb_filenames]

    s_widget = IntSlider(min=0, max=len(pdbstrings)-1, description="Structure", continuous_update=False)
    return interact(view_py3Dmol, i=s_widget)


def view_pose(pose=None, hbonds=True, residue_selector=None, label=False):
    """View a PyRosetta ResidueSelector in a PyRosetta Pose object in py3Dmol within a Jupyter notebook.
    Optionally also label the residues in the PyRosetta ResidueSelector or show all hydrogen bonds.
    The user must have already initialized PyRosetta providing .params files for any ligands/non-canonical residues
    in the pose.
    """
    if not isinstance(pose, pyrosetta.rosetta.core.pose.Pose):
        raise ValueError("Input pose should be an instance of pyrosetta.rosetta.core.pose.Pose")
    
    viewer = py3Dmol.view(1200, 800)
    viewer.addModels(pyrosetta.distributed.io.to_pdbstring(pose), "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}, "stick": {"colorscheme": "blackCarbon", "radius": 0.025}})

    if hbonds:
        hbond_set = pose.get_hbonds()
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_xyz = pose.residue(r.don_res()).xyz(r.don_hatm())
                    acc_xyz = pose.residue(r.acc_res()).xyz(r.acc_atm())
                    viewer.addLine({"dashed": True, "color": "black",
                                    "start": {"x": don_xyz[0], "y": don_xyz[1], "z": don_xyz[2]},
                                    "end": {"x": acc_xyz[0], "y": acc_xyz[1], "z": acc_xyz[2]}})

    if residue_selector:
        if not isinstance(residue_selector, pyrosetta.rosetta.core.select.residue_selector.ResidueSelector):
            raise ValueError(
                "Input residue_selector should be an instance of pyrosetta.rosetta.core.select.residue_selector.ResidueSelector"
            )   
        residue_list = list(pyrosetta.rosetta.core.select.get_residues_from_subset(residue_selector.apply(pose))) 
        viewer.setStyle({"resi": residue_list},
                        {"cartoon": {"color": "spectrum"}, "stick": {"colorscheme": "whiteCarbon", "radius": 0.25}})  
        if label:
            viewer.addResLabels({"resi": residue_list}, {"fontSize": 12, "showBackground": False, "fontColor": "black"})
    
    viewer.zoomTo()
    return viewer.show()
