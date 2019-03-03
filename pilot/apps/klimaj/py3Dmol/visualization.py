import os
import py3Dmol
import pyrosetta
import pyrosetta.distributed.io
import pyrosetta.distributed.packed_pose

from ipywidgets import interact as _interact
from ipywidgets import IntSlider as _IntSlider
from IPython.core.display import display as _display
from IPython.core.display import HTML as _HTML

__author__ = "Jason C. Klima"
__email__  = "klimaj@uw.edu"


_display(_HTML("<style>.container { width:100% !important; }</style>"))

def view_poses(poses=None, window_size=(1200, 800)):
    """View a list of PyRosetta Pose or PackedPose objects in py3Dmol with the ipywidgets
    interactive slider within a Jupyter notebook.
    Optionally change the displayed window size. The user must have already initialized PyRosetta
    providing .params files for any ligands and non-canonical residues in the poses.
    """

    def view_py3Dmol(i=0):
        viewer = py3Dmol.view(*window_size)
        viewer.addModels(pdbstrings[i], "pdb")
        viewer.setStyle(
            {"cartoon": {"color": "spectrum"}, "stick": {"radius": 0.25}}
        )
        viewer.setStyle(
            {"hetflag": True},
            {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}}
        )
        viewer.zoomTo()

        return viewer.show()


    if not isinstance(poses, list):
        raise ValueError("Input should be a list of poses and/or packed poses.")
    else:
        for p in poses:
            if isinstance(p, pyrosetta.rosetta.core.pose.Pose):
                pass
            elif isinstance(p, pyrosetta.distributed.packed_pose.core.PackedPose):
                pass
            else:
                raise ValueError("Input poses should be instances of pyrosetta.rosetta.core.pose.Pose \
                    or pyrosetta.distributed.packed_pose.core.PackedPose")

    if not isinstance(window_size, tuple):
       raise ValueError("Input window size should be an instance of tuple.")
    assert len(window_size) == 2, "Input window size must be a tuple of length 2."

    pdbstrings = [
        pyrosetta.distributed.io.to_pdbstring(p) for p in poses
    ]
    s_widget = _IntSlider(
        min=0, max=len(pdbstrings)-1, description="Decoys", continuous_update=False
    )
    
    return _interact(view_py3Dmol, i=s_widget)


def view_pdbs(pdb_filenames=None, window_size=(1200, 800)):
    """View a list of .pdb files in py3Dmol with the ipywidgets interactive slider within a Jupyter
    notebook.
    Optionally change the displayed window size. The user must have already initialized PyRosetta
    providing .params files for any ligands and non-canonical residues in the .pdb files.
    """

    def view_py3Dmol(i=0):
        viewer = py3Dmol.view(*window_size)
        viewer.addModels(pdbstrings[i], "pdb")
        viewer.setStyle(
            {"cartoon": {"color": "spectrum"}, "stick": {"radius": 0.25}}
        )
        viewer.setStyle(
            {"hetflag": True},
            {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}}
        )
        viewer.zoomTo()

        return viewer.show()


    if not isinstance(pdb_filenames, list):
        raise ValueError("Input should be a list of valid .pdb file paths.")
    else:
        for pdb_filename in pdb_filenames:
            if not os.path.exists(pdb_filename):
                raise ValueError("Input .pdb files should be valid path strings.")

    if not isinstance(window_size, tuple):
        raise ValueError("Input window size should be an instance of tuple.")
    assert len(window_size) == 2, "Input window size must be a tuple of length 2."

    pdbstrings = [
        pyrosetta.distributed.io.to_pdbstring(
            pyrosetta.distributed.io.pose_from_file(p)
        ) for p in pdb_filenames
    ]
    s_widget = _IntSlider(
        min=0, max=len(pdbstrings)-1, description="Decoys", continuous_update=False
    )
    
    return _interact(view_py3Dmol, i=s_widget)


def view_pose(pose=None, residue_selector=None, labels=True,
              hbonds=True, disulfides=True, window_size=(1200, 800)):
    """View a PyRosetta Pose or PackedPose object in py3Dmol within a Jupyter notebook.
    Also optionally view a PyRosetta ResidueSelector with or without labeling the residues in the
    PyRosetta ResidueSelector using PDB numbering, optionally show all hydrogen bonds according to
    PyRosetta, optionally show all disulfide bonds according to PyRosetta, and optionally change the
    displayed window size. The user must have already initialized PyRosetta providing .params files
    for any ligands and non-canonical residues in the pose.
    """

    if isinstance(pose, pyrosetta.rosetta.core.pose.Pose):
        pass
    elif isinstance(pose, pyrosetta.distributed.packed_pose.core.PackedPose):
        pose = pyrosetta.distributed.packed_pose.to_pose(pose)
    else:
        raise ValueError("Input pose should be an instance of pyrosetta.rosetta.core.pose.Pose \
            or pyrosetta.distributed.packed_pose.core.PackedPose")

    if not isinstance(window_size, tuple):
        raise ValueError("Input window size should be an instance of tuple.")
    assert len(window_size) == 2, "Input window size must be a tuple of length 2."

    viewer = py3Dmol.view(*window_size)
    viewer.addModels(pyrosetta.distributed.io.to_pdbstring(pose), "pdb")

    if residue_selector:

        if not isinstance(residue_selector,
                          pyrosetta.rosetta.core.select.residue_selector.ResidueSelector):
            raise ValueError("Input residue_selector should be an instance of \
                pyrosetta.rosetta.core.select.residue_selector.ResidueSelector"
            )

        major_radius = 0.5
        minor_radius = 0.05

        viewer.setStyle(
            {"cartoon": {"color": "spectrum"},
            "stick": {"colorscheme": "blackCarbon", "radius": minor_radius}}
        )

        pdb_numbering = list(zip(
                    pyrosetta.rosetta.core.pose.full_model_info.get_res_num_from_pdb_info(pose),
                    pyrosetta.rosetta.core.pose.full_model_info.get_chains_from_pdb_info(pose)
        ))
        residues_from_subset = list(
            pyrosetta.rosetta.core.select.get_residues_from_subset(residue_selector.apply(pose))
        )
        residue_chain_tuples = [pdb_numbering[i-1] for i in residues_from_subset]
        
        for resi, chain in residue_chain_tuples:
            viewer.setStyle(
                {"resi": resi, "chain": chain}, {"cartoon": {"color": "spectrum"},
                "stick": {"colorscheme": "whiteCarbon", "radius": major_radius}}
            )
        
        if labels:
            for resi, chain in residue_chain_tuples:
                viewer.addResLabels(
                    {"resi": resi, "chain": chain},
                    {"fontSize": 12, "showBackground": False, "fontColor": "black"}
                )
                
    else:

        minor_radius = 0.25
        viewer.setStyle(
            {"cartoon": {"color": "spectrum"},
            "stick": {"colorscheme": "lightgreyCarbon", "radius": minor_radius}}
        )

    if hbonds:

        hbond_set = pose.get_hbonds()
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_xyz = pose.residue(r.don_res()).xyz(r.don_hatm())
                    acc_xyz = pose.residue(r.acc_res()).xyz(r.acc_atm())
                    viewer.addLine(
                        {
                            "dashed": True, "color": "black",
                            "start": {"x": don_xyz[0], "y": don_xyz[1], "z": don_xyz[2]},
                            "end": {"x": acc_xyz[0], "y": acc_xyz[1], "z": acc_xyz[2]}
                        }
                    )

    if disulfides:

        cys_res = []
        for i, aa in enumerate(pose.sequence(), start=1):
            if aa == "C":
                cys_res.append(i)
        for i in cys_res:
            for j in cys_res:
                if pyrosetta.rosetta.core.conformation.is_disulfide_bond(pose.conformation(), i, j):
                    i_xyz = pose.xyz(
                        pyrosetta.rosetta.core.id.AtomID(pose.residue(i).atom_index("SG"), i)
                    )
                    j_xyz = pose.xyz(
                        pyrosetta.rosetta.core.id.AtomID(pose.residue(j).atom_index("SG"), j)
                    )
                    viewer.addCylinder(
                        {
                            "radius": minor_radius, "color": "gold", "fromCap": True, "toCap": True,
                            "start": {"x": i_xyz[0], "y": i_xyz[1], "z": i_xyz[2]},
                            "end": {"x": j_xyz[0], "y": j_xyz[1], "z": j_xyz[2]}
                        }
                    ) 
    
    viewer.zoomTo()

    return viewer.show()
