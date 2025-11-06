from pyrosetta import *
import numpy as np
from bootcamp_mover import BootCampMover

class BootCampMoverCreator(rosetta.protocols.moves.MoverCreator):

    _instances = list()

    def __init__(self):
        rosetta.protocols.moves.MoverCreator.__init__(self)

    def create_mover(self):
        mover = BootCampMover()
        self._instances.append(mover)
        return mover

    def keyname(self):
        return BootCampMover.mover_name()

    def provide_xml_schema(self, xsd):
        print("Creator provide_xml_schema is called. ")
        BootCampMover.provide_xml_schema(xsd)


#global variable
_py_mover_creators_ = []

def register():
    factory = rosetta.protocols.moves.MoverFactory.get_instance()
    creator = BootCampMoverCreator()
    factory.factory_register(creator)

    _py_mover_creators_.append(creator)
