"""
--------------------------------------------------------------------------------
 Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 This file is part of FLEUR and available as free software under the conditions
 of the MIT license as expressed in the LICENSE file in more detail.
--------------------------------------------------------------------------------

DFPT Phonon Data Reader and JSON Export Utilities
=================================================

This module provides a collection of tools for reading, processing, and 
serializing phonon data obtained from Density Functional Perturbation Theory 
(DFPT) calculations—specifically, the `dynMatq*` files produced by FLEUR.

The module covers the needed functions to prepare phonon bandstructure
and eigenvector data for visualization in external tools (e.g. web-based 
phonon viewers, Python dashboards, or custom GUI applications).

Main Capabilities
-----------------
1. **Reading DFPT output files**
   - `read_dynMats`           → read 3N×3N dynamical matrices for all q-points  
   - `read_eigenvectors`      → extract phonon eigenvectors (complex)  
   - `read_eigenfrequencies`  → extract phonon frequencies and q-vectors  

2. **Coordinate and format conversions**
   - `cast_cartesian`         → convert fractional q-points to Cartesian space  
   - `calc_q_distance`        → compute cumulative distances along a q-path  
   - `recast_eigenvecs_dim`   → reshape eigenvectors for visualization  
   - `recast_eigenmodes`      → flatten eigenmodes into (3N) representation  
   - `get_filename`           → construct correctly zero-padded dynMat filenames  

3. **Packaging data for visualization**
   - `generate_dict`          → build a structured dictionary containing all
                                phonon information (frequencies, eigenvectors,
                                lattice, atomic positions, q-path metadata, etc.)
   - `save_json`              → export the dictionary to JSON using a custom
                                encoder that supports NumPy arrays and complex
                                numbers

4. **JSON Serialization Support**
   - `NumpyEncoder`           → safely serializes NumPy arrays, scalar types,
                                and complex numbers into JSON-compatible formats
"""


import numpy as np
import json




def generate_dict(
    distances,
    natoms,
    eigenvectors,
    system_name,
    eigenvalues,
    repetitions,
    chemical_symbols,
    qpoints,
    atom_numbers,
    lattice,
    #atomic_numbers,
    highsym_qpts,
    atom_pos_red,
    atom_pos_car,
    formula,
    atom_types
):
    """
   Generate a dictionary containing all relevant data for a DFPT Phonon calculation. This dict can be exported for visualization

   Parameters
   ----------
   distances : list or array-like 
       Distances between the q-points 
   natoms : int
       Number of atoms in the system.
   vectors : array-like (Nq,3*natoms,natoms,3,2)
       Eigenvectors of the system 
   system_name : str
       Name or label for the system.
   eigenvalues : list or array-like (Nq,3*natoms)
       Eigenvalues of the dynamical matrix at each q-point 
   repetitions : list or array-like (3)
       Supercell repetitions along each lattice direction.
   chemical_symbols : list of str
       Chemical symbols for each atom in the system.
   qpoints : list or array-like (Nq,3)
       List of q-points used in the phonon calculations.
   atom_numbers : list of int
       Atomic numbers for each atom.
   lattice : array-like (3,3)
       Bravais Matrix of the system. 
   highsym_qpts : list or array-like
       Labels for q-points along the path.
   atom_pos_red : array-like
       Atomic positions in reduced (fractional) coordinates.
   atom_pos_car : array-like
       Atomic positions in Cartesian coordinates.
   formula : str
       Chemical formula of the system.
   atom_types : list of str or int
       Symbol of the elements in the system
       --> this seems redundant (chemical_symbols) but 
           needed for the visualisation. 

   Returns
   -------
   dict
       A dictionary containing all the input data
   """
    data = {
        "distances": distances,
        "natoms": natoms,
        "vectors": eigenvectors,
        "name": system_name,
        "eigenvalues": eigenvalues,
        "repetitions": repetitions,
        "chemical_symbols": chemical_symbols,
        "qpoints": qpoints,
        "atom_numbers": atom_numbers,
        "lattice": lattice,
        #"atomic_numbers": atomic_numbers,
        "highsym_qpts": highsym_qpts,
        "atom_pos_red": atom_pos_red,
        "atom_pos_car": atom_pos_car,
        "formula": formula,
        "atom_types": atom_types
    }
    return data


def save_json(filename, data):
    with open(filename, "w") as f:
        json.dump(data, f, indent=4, cls=NumpyEncoder)
    

def read_dynMats(path,number_atoms,number_q,prefix=''):
    """
    Read dynamical matrices from a FLEUR DFPT calculation. 

    This function loads the 3N × 3N dynamical matrix (where N is the number
    of atoms) for each q-point. The matrices are read from plain-text files
    produced by a FLEUR DFPT calculations (`dyn*` files).
    Each file is assumed to follow the standard format where real and
    imaginary parts of the matrix elements alternate in the columns.

    Parameters
    ----------
    path : str
        Path to the folder that contains the dynMats. 
        
    number_atoms : int
        Number of atoms in the unit cell (N).
    number_q : int
        Number of q-points 
    prefix : str, optional
        String to prepend to each generated filename. 
        Needed if band interpolation is used (e.g. band_)

    Returns
    -------
    dyn : ndarray, shape (number_q, 3N, 3N)
        A NumPy array containing the complex dynamical matrices for all q-points.
        - `dyn[q]` corresponds to the dynamical matrix at q-point index `q`
        - Each matrix element is complex: real part + i * imaginary part
    """
    dyn = np.zeros((number_q,3*number_atoms,3*number_atoms),dtype=complex)
    for q in range(1,number_q+1):
        file_name = get_filename(path,q,prefix)
        with open(file_name,"r") as file:
            lines = file.readlines()
            counter_row = 0# runs --> 
            counter_col = 0 # runs vertically 
            for line in lines[4:]:
                cont = line.split()
                tmp_r = cont[0::2]
                tmp_c = cont[1::2]
                for i in range(3):
                    dyn[q-1][counter_col][counter_row] = float(tmp_r[i]) + 1j* float(tmp_c[i])
                    counter_row+=1
                    if counter_row == 3*number_atoms:
                        counter_row = 0
                        counter_col +=1

                if counter_col == 3*number_atoms:
                    break  # full dynmat read

    return dyn


def read_eigenvectors(path,number_q,number_atoms,prefix=''):
    """
    Read phonon eigenvectors from DFPT output files.

    This function loads the complex phonon eigenvectors for each q-point from
    text files dynMat*. 

    The function parses the file past the dynamical matrix and additional
    header lines, then reads the eigenvector blocks until the
    "Eigenfrequencies" marker or until the full eigenvector matrix is filled.

    Parameters
    ----------
    path : str
        path to the folder that contains the dynMatq* files.
    number_q : int
        Number of q-points.
    number_atoms : int
        Number of atoms in the system (N).
    prefix : str, optional
        String to prepend to each generated filename. 
        Needed if band interpolation is used (e.g. band_)

    Returns
    -------
    ndarray
        A complex NumPy array containing the eigenvectors with shape:
        `(number_q, 3N, N,3,2)` 
        Complex and realpart are splitted in the last index for the .json format

    """
    eigvec = np.zeros((number_q,3*number_atoms,3*number_atoms),dtype=complex)
    for q in range(1,number_q+1):
        file_name = get_filename(path,q,prefix)
        with open(file_name,"r") as file:
            lines = file.readlines()
            counter_row = 0# runs --> 
            counter_col = 0 # runs vertically 
            
            start_line = 4 + (3*number_atoms*3*number_atoms)/3 + number_atoms
            for line in lines[int(start_line):]:
                cont = line.split()
                if line.strip() == "":
                    continue  # Skip empty lines
                try:
                    if cont[0] == "Eigenfrequencies":
                        break 
                except IndexError:
                    continue 
                tmp_r = cont[0::2]
                tmp_c = cont[1::2]
                for i in range(3):
                    eigvec[q-1][counter_col][counter_row] = float(tmp_r[i]) + 1j* float(tmp_c[i])
                    counter_row+=1
                    if counter_row == 3*number_atoms:
                        counter_row = 0
                        counter_col +=1

                if counter_col == 3*number_atoms:
                    break  # full dynmat read

    return recast_eigenvecs_dim(eigvec, number_q, number_atoms)



def read_eigenfrequencies(path,number_q,number_atoms,prefix=''):
    """
   Read phonon eigenfrequencies from FLEUR DFPT dynMat files.

   For each q-point, it extracts:

   1. The q-vector (fractional coordinates)
   2. The phonon eigenfrequencies:
      - real and imaginary components

   The function returns the q-path and an eigenfrequencies array. 

   Parameters
   ----------
   path : str
       path to the dynMat files. 
   number_q : int
       Number of q-points 
   number_atoms : int
       Number of atoms N in the system.
   prefix : str, optional
        String to prepend to each generated filename. 
        Needed if band interpolation is used (e.g. band_)

   Returns
   -------
   q_sample : ndarray, shape (number_q, 3)
       Fractional coordinates of each q-point along the phonon path.
   eigenmodes : ndarray, , shape (number_q, 3N)
       Complex eigenfrequency arry. 
       
   """
    q_sample=np.zeros((number_q,3)) #Store all q vectors along the path
    eigenfrequencies=np.zeros((number_q,3,number_atoms)) #store all eigenmodes in a q_N times 3*Atoms array

    for i in range(1,number_q+1):
        file_name = get_filename(path,i,prefix)
        with open(file_name,'r') as file:
            lines=file.readlines()
            for j in range(0,len(lines)):    
                line=lines[j]
                if "q =" in line:
                    q = [ float(x) for x in line.split()[2:] ]
                    q_sample[i-1,:]=np.array(q)
                if "Eigenfrequencies in 1/cm" in line:
                    l_next=lines[j+1]
                    eigenfreq=[ float(x) for x in l_next.split()[2:-1 ] ] #obtain the  eigenfreq 1/cm with entries Re(a1), Im(a1) , Re(a2))
                    eigenfrequencies[i-1,:,0]=np.array(eigenfreq[0::2]) #slicing such that we only get the real part of the eigenmodes 
                    eigenfrequencies[i-1,:,0]+=np.array(eigenfreq[1::2]) #slicing to get the imaginary part 
                    if number_atoms>1:
                        for t_atoms in range(1,number_atoms):
                            l_type=lines[j+1+t_atoms]
                            eigenfreq=[ float(x) for x in l_type.split()[2:-1 ] ] #obtain the  eigenfreq 1/cm with entries Re(a1), Im(a1) , Re(a2))
                            eigenfrequencies[i-1,:,t_atoms]=np.array(eigenfreq[0::2]) #slicing such that we only get the real part of the eigenmodes
    
        
    return q_sample,recast_eigenmodes(eigenfrequencies,number_q,number_atoms)


def calc_q_distance(arr,bmat,number_q):
    """
    Compute cumulative distances along a q-point path in reciprocal space.
 
    This function converts fractional q-points to Cartesian coordinates using
    the reciprocal lattice matrix `bmat`, then computes the cumulative
    distance along the q-path. The result is stored in distance
 
    Parameters
    ----------
    arr : ndarray, shape (number_q, 3)
        List of q-points in reduced (fractional) coordinates.
    distance : ndarray, shape (number_q)
        Output array that will be filled with cumulative distances.
        - `distance[0]` is incremented by the distance to the BZ-zone center
        - `distance[q]` contains the cumulative arc length of the path.
    bmat : ndarray, shape (3, 3)
        Reciprocal lattice matrix.
    number_q : int
        Number of q-points in the path.
 
    Returns
    -------
    distance
   
    """    
    q_cart = cast_cartesian(arr, bmat)
    
    distance = np.zeros(number_q)

    distance[0] += np.sqrt(q_cart[0,0]**2+q_cart[0,1]**2+q_cart[0,2]**2)
    
    for q in range(1,number_q):
        dist = q_cart[q,:] - q_cart[q-1,:]
        cart = np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
        distance[q] = distance[q-1] + cart  

    return distance
    




"""
Helper functions
"""

class NumpyEncoder(json.JSONEncoder):
    """
    JSON encoder that safely serializes NumPy data types.
   
    This custom JSON encoder extends Python's built-in ``json.JSONEncoder`` to
    allow serialization of commonly used NumPy objects that are not JSON
    serializable by default. It handles:
   
    - ``np.ndarray`` → converted to nested Python lists
    - ``np.integer`` → converted to Python ``int``
    - ``np.floating`` → converted to Python ``float``
    - complex numbers → encoded as ``[real, imag]`` pairs
   
    """
    def default(self, obj):
        # NumPy arrays → lists
        if isinstance(obj, np.ndarray):
            return obj.tolist()

        # NumPy integer / float scalars
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)

        # Complex numbers → [real, imag]
        if isinstance(obj, complex):
            return [obj.real, obj.imag]

        return super().default(obj)
 



def recast_eigenvecs_dim(arr,number_q,number_atoms):
    
    """
    Recast eigenvector matrix into a structured multidimensional format.

    This function takes the raw DFPT eigenvector array of shape
    ``(number_q, 3N, 3N)``—where N is the number of atoms—and reorganizes it
    into a structure that can be used for visualization. 

        (number_q, 3N, N, 3, 2)
        
    """
    
    # we save eigenvectors as [nq,dir,vec]

    recast = np.zeros((number_q,3*number_atoms,number_atoms,3,2))


    for iSuperDir in range(3*number_atoms):
        for iVec in range(3*number_atoms):
            element = np.copy(arr[:,iSuperDir,iVec])
            
            iAtom = iSuperDir // 3     # integer division
            iDir  = iSuperDir % 3      # remainder
        
            
                
            recast[:,iVec,iAtom,iDir,0] = np.real(element)
            recast[:,iVec,iAtom,iDir,1] = np.imag(element)
        
    return recast
        
def recast_eigenmodes(arr,number_q,number_atoms):
    """
    Recast eigenmode data into a flattened (3N) phonon-branch representation.

    This function takes phonon eigenmodes stored in the format
    ``(number_q, 3, number_atoms)`` and 
    rearranges them into a flattened array of shape:

        (number_q, 3N)

    where N is the number of atoms, and the 3 Cartesian components of each atom
    are placed contiguously:

        [u_x(atom0), u_y(atom0), u_z(atom0), u_x(atom1), ..., u_z(atomN−1)]
    """
    
    # we save eigenvecotrs as [nq,dir,vec]

    recast = np.zeros((number_q,3*number_atoms))


    for iAtom in range(number_atoms):
            element = np.copy(arr[:,:,iAtom])
            
            #iAtom = iSuperDir // 3     # integer division
            #iDir  = iSuperDir % 3      # remainder
        
            start = 3 * iAtom
            end   = 3 * (iAtom + 1)
                
            recast[:,start:end] = element 
        
        
    return recast


def cast_cartesian(arr,bmat):
    """
    Convert reduced (fractional) coordinates to Cartesian coordinates.

    Given an array of q-points (or atomic positions) expressed in reduced
    coordinates and a reciprocal lattice matrix ``bmat``, this function computes
    their Cartesian representation:

        q_cart = bmat · q_frac
    """
    y = np.zeros(np.shape(arr))
    for i in range(len(arr[:,0])):
        y[i,:]=np.matmul(bmat,np.transpose(arr[i,:]))
        
    return y




def get_filename(f,q,prefix=''):
    """
    Build the filename for the dynmMat.

        dynMatq=0001
        dynMatq=0012
        dynMatq=0123
        dynMatq=1234

    This function formats the q-point index with the correct zero-padding
    """
    if q<10:
        filename = f'{f}/{prefix}dynMatq=000{q}'
    elif q>=10 and q<100:
        filename = f'{f}/{prefix}dynMatq=00{q}'
    elif q>=100 and q<1000:
        filename = f'{f}/{prefix}dynMatq=0{q}'
    else:
        filename = f'{f}/{prefix}dynMatq={q}'

    return filename


