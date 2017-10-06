# -*- coding: utf-8 -*-
"""
Module to work on Ane's problem of efficiency in Tractography streamlines
"""
import nibabel as nib
import numpy as np
import os
from os.path import join as opj
import subprocess

CWD = os.getcwd()
"""
1: Get Corpus Callosum mask and extract the upper center part.
ie: fix [y] and [z] (upper points in median sagittal plane)
This will be the ground truth mask
"""


def create_corpus_callosum():

    # JHU DTI-based white-matter atlases
    JHU = '/usr/share/fsl/5.0/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz'

    JHU_img = nib.load(JHU)
    atlas_data = JHU_img.get_data()
    corpus_callosum_data = np.zeros((atlas_data.shape))

    # Genu of Corpus Callosum
    corpus_callosum_data[np.where(atlas_data == 3)] = 1
    # Body of Corpus Callosum
    corpus_callosum_data[np.where(atlas_data == 4)] = 1
    # Splenium of Corpus Callosum
    corpus_callosum_data[np.where(atlas_data == 5)] = 1

    corpus_callosum_img = nib.Nifti1Image(corpus_callosum_data,
                                          affine=JHU_img.affine)

    nib.save(corpus_callosum_img, opj(CWD,
                                      'data',
                                      'corpus_callosum_1mm.nii.gz'))


def create_corpus_callosum_plane():
    # mid-sagital plane
    from nilearn.plotting import plot_glass_brain
    # JHU DTI-based white-matter atlases
    corpus_callosum = opj(CWD, 'data', 'corpus_callosum_1mm.nii.gz')

    corpus_callosum_img = nib.load(corpus_callosum)
    corpus_callosum_data = corpus_callosum_img.get_data()
    corpus_callosum_med_sag_plane = np.zeros((corpus_callosum_data.shape))

    # mid plane-1 to fit with MNI(x)=0
    mid_sag_plane = (corpus_callosum_data.shape[0]//2)-1

    corpus_callosum_med_sag_plane[mid_sag_plane,:,:] = corpus_callosum_data[mid_sag_plane,:,:]

    corpus_callosum_med_sag_plane_img = nib.Nifti1Image(corpus_callosum_med_sag_plane,
                                                        affine=corpus_callosum_img.affine)

    nib.save(corpus_callosum_med_sag_plane_img, opj(CWD,
                                                    'data',
                                                    'corpus_callosum_med_sag_plane_1mm.nii.gz'))    



def create_skeleton_atlas():
    # Just area/surface that actually can be passed trough
    # TODO: This is more difficult than I though
    
    # For each point of the MNI atlas mask, calculate the min distance to 
    # CC_med_sag_plane. Create a line of points between the points and check 
    # for points that in CC x,y,z boundaries, fall out of CC.



from scipy.spatial.distance import cdist

def closest_node(node, nodes):
    return nodes[cdist([node], nodes).argmin()]

"""
2: Transform the CC mask to subjects DWI space
"""


def execute(cmd):
    popen = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def transform_mask_to_subject_space(mask_path=MASK_PATH,
                                    t1_atlas_path=MNI_1MM_PATH,
                                    dwi_subject_path)

    import tempfile
    
    # review this naming:
    mask_dwi_path = dwi_subject_path[:-7]
    
    if os.path.exists(mask_dwi_path):
        return
    
    omat = tempfile.mkstemp()

    if not os.path.exists(omat):
        command = ['flirt',
                   '-in',
                   t1_atlas_path,
                   '-ref',
                   dwi_subject,
                   '-omat',
                   omat[1],
                   ]
        for output in execute(command):
            print(output)

    command = ['flirt',
               '-in',
               mask_path,
               '-ref',
               dwi_subject_path,
               '-out',
               mask_dwi_path,
               '-init',
               omat[1],
               '-applyxfm', '-interp', 'nearestneighbour',
               ]
    for output in execute(command):
        print(output)







"""
3: Calculate the most efficient path between 2 given points in the space 
crossing CC mask
"""







