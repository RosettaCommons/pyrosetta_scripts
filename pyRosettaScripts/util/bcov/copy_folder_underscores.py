#!/usr/bin/python


import shutil
import sys
import os
import subprocess

src = sys.argv[1]
dest = sys.argv[2]




def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return the_stuff[0] + the_stuff[1]



new_name = src.replace("/", "_")

dest = os.path.join(dest, new_name)

print cmd("rsync -avL %s %s"%(src, dest))


