import os
import pyrosetta
import pyrosetta.toolbox
import pyrosetta.distributed.dask

from dask import delayed
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, fire_and_forget


__author__ = "Jason C. Klima"
__email__  = "klimaj@uw.edu"


def my_protocol(start_pose):
    scorefxn = pyrosetta.create_score_function("ref2015")
    protocol = pyrosetta.rosetta.protocols.moves.NullMover()
    jd = pyrosetta.toolbox.py_jobdistributor.PyJobDistributor(pdb_name=start_pose.pdb_info().name().split(".pdb")[0],
                                                              nstruct=10,
                                                              scorefxn=scorefxn,
                                                              compress=True)
    jd.native_pose = start_pose
    jd.json_format = True
    jd.additional_decoy_info = "Hello from SLURM_JOB_ID " + os.environ["SLURM_JOB_ID"]
    while not jd.job_complete:
        pose = start_pose.clone()
        protocol.apply(pose)
        jd.output_decoy(pose)
    return start_pose


def main(pdbfile, flags):
    pyrosetta.distributed.dask.init_notebook(flags)
    scratch_dir = os.path.join("/net/scratch", os.environ["USER"])
    cluster = SLURMCluster(cores=1,
                           processes=1,
                           memory="4GB",
                           queue="short",
                           local_directory=scratch_dir,
                           job_extra=["-o {}".format(os.path.join(scratch_dir, "slurm-%j.out"))],
                           extra=pyrosetta.distributed.dask.worker_extra(init_flags=flags))
    cluster.scale(1)
    client = Client(cluster)
    start_pose = pyrosetta.io.pose_from_file(pdbfile)
    fire_and_forget(delayed(my_protocol)(start_pose).compute())


if __name__ == "__main__":
    pdbfile = "crambin.pdb"
    flags = """
    -ignore_unrecognized_res 1
    -ex4
    """
    main(pdbfile=pdbfile, flags=flags)
