# -*- coding: utf-8 -*-
"""
Module to work on Ane's problem of efficiency in Tractography streamlines
"""
import nibabel as nib
import numpy as np
import os
from os.path import join as opj

from scipy.spatial.distance import cdist

from efficiency.utils import (execute,
                              closest_node,
                              bresenhamline,
                              )

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
    
    # Full CC
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

# discuss if we need this or not:

#def create_not_corpus_callosum_vol():
#    
#    # Full CC
#    corpus_callosum = opj(CWD, 'data', 'corpus_callosum_1mm.nii.gz')
#
#    corpus_callosum_img = nib.load(corpus_callosum)
#    corpus_callosum_data = corpus_callosum_img.get_data()
#
#    MNI_brain = '/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain_mask.nii.gz'  
#    MNI_brain_data = nib.load(MNI_brain).get_data()
#    
#    not_corpus_callosum = MNI_brain_data - 
#
#    corpus_callosum_med_sag_plane[mid_sag_plane,:,:] = corpus_callosum_data[mid_sag_plane,:,:]
#
#    corpus_callosum_med_sag_plane_img = nib.Nifti1Image(corpus_callosum_med_sag_plane,
#                                                        affine=corpus_callosum_img.affine)
#
#    nib.save(corpus_callosum_med_sag_plane_img, opj(CWD,
#                                                    'data',
#                                                    'corpus_callosum_med_sag_plane_1mm.nii.gz'))   

"""
2: Transform the CC mask to subjects DWI space
"""


def transform_mask_to_subject_space(mask_path=MASK_PATH,
                                    t1_atlas_path=MNI_1MM_PATH,
                                    dwi_subject_path):

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

def find_optimal_cc_crossing(point_1, point_2, area):
    """
    Function to calculate the optimal point of cut between 2 points given the
    constraint of having to cross an area of points.

    Returns min_distance and optimal_point
    
    Example:
        
    # JHU DTI-based white-matter atlases
    corpus_callosum_med_sag =  opj(CWD,
                                   'data',
                                   'corpus_callosum_med_sag_plane_1mm.nii.gz')
        
    corpus_callosum_med_sag_img = nib.load(corpus_callosum_med_sag)
    area = corpus_callosum_med_sag_img.get_data()
    
    point_1 = np.array([24, 45, 45])
    point_2 =  np.array([120, 90, 90])
    
    >>> find_optimal_cc_crossing(point_1, point_2, area)
    >>> (116.58400271916257, array([90, 84, 81]))
    """

    x, y, z = np.where(area == 1)
    area_points = np.array([[x, y, z] for x, y, z in zip(x, y, z)])

    # Calculating all options (it is expected not the be much area),
    # and find the minimum.
    distances = [cdist([area_point], [point_1]) +
                 cdist([area_point], [point_2])
                 for area_point in area_points]

    min_distance = float(min(distances))
    optimal_point = area_points[argmin(distances)]

    return min_distance, optimal_point

