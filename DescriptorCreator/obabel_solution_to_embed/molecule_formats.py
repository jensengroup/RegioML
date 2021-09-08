import subprocess

def shell(cmd, shell=False):

    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, err = p.communicate()
    return output


def convert_sdf_to_xyz(sdffile, xyzfile):
    
    shell(f'obabel -isdf {sdffile} -oxyz -xf > {xyzfile}', shell=True)

    return

def convert_2Dsdf_to_3Dxyz(sdffile, xyzfile):
    
    shell(f'obabel -isdf {sdffile} -oxyz --gen3D > {xyzfile}', shell=True)

    return
