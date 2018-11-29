import sys
from os.path import join as pjoin
import os


class Slurm(object):
    """docstring for Slurm."""

    ssh_command = "ssh -t {user}@{server} '#!/bin/bash\nsbatch {script}'"

    scp_command = 'scp {files} {user}@{server}:{folder}'

    back_base = '{user}@{server}:{folder}/{file}'

    scp_back = 'scp -r {files} /tmp/'

    slurm_script = """#!/bin/bash
#SBATCH -D {folder}
#SBATCH -J {job_name}
#SBATCH -o {folder}/{job_name}.out
#SBATCH -e {folder}/{job_name}.err
#SBATCH -A {project}
#SBATCH -t {time}
#SBATCH -n {threads}
#SBATCH -p core

{commands}
"""

    def run_remote(self, job_name, command, files, threads = 20, folder = "/home/moritz/temp/", time = '10-00:00:00', project = "snic2017-1-616", email = "murumbii@google.com"):
        script = self.slurm_script.format(
            folder = folder,
            job_name = job_name,
            project = project,
            time = time,
            threads = threads,
            commands = command)

        with open(pjoin('/tmp/', job_name + ".sh"), "w") as  handle:
            handle.writelines(script)

        files = [pjoin('/tmp/', job_name + ".sh")] + files
        scp = self.scp_command.format(files = " ".join(files), user = "moritz", server = "rackham.uppmax.uu.se", folder = folder )
        ssh = self.ssh_command.format(user = "moritz", server = "rackham.uppmax.uu.se", script = pjoin(folder, job_name + ".sh"))

        os.system(scp)
        os.system(ssh)
        self.clean_up(files)

    def retrieve_data(self, files, folder):
        file_full = [self.back_base.format(file = f, user = "moritz", server = "rackham.uppmax.uu.se", folder = folder) for f in files]

        scp = self.scp_back.format(files = " ".join(file_full))

        os.system(scp)

    def clean_up(self, files):
        for f in files:
            os.remove(f)
