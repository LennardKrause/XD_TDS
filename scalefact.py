import os
import re
import numpy as np
import pandas as pd

def get_cell_and_wave(fpath):
    """
    Extracts unit cell parameters and wavelength from a file.

    Parameters:
        fpath (str): Path to the input file containing 'CELL' and 'WAVE' records.

    Returns:
        tuple:
            cell (numpy.ndarray): Array of six floats representing the unit cell parameters (a, b, c, alpha, beta, gamma).
            wave (float): The wavelength value extracted from the file.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        IndexError: If 'CELL' or 'WAVE' records are not found in the file.
        ValueError: If the extracted values cannot be converted to float.
    """
    with open(fpath, 'r') as f:
        content = f.read()
    cell = re.findall(r'CELL\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', content)[0]
    cell = np.array(list(map(float, cell)))
    wave = re.findall(r'WAVE\s+(\d+\.\d+)', content)[0]
    wave = float(wave)
    return cell, wave

def get_data(fpath):
    """
    Reads a data file and parses its contents into a pandas DataFrame.

    Parameters:
        fpath (str): Path to the data file to be read.

    Returns:
        tuple:
            head (str): The first line (header) of the file.
            hkl (pandas.DataFrame): DataFrame containing the parsed data with columns named according to the file format.
            ndat (int): Number of data columns as determined from the header.

    Raises:
        ValueError: If the number of data columns in the header is not 6 or 7.
    """
    with open(fpath, 'r') as f:
        head = f.readline()
    ndat = int(head.split()[-1].strip())
    if ndat == 7:
        names = ['h', 'k', 'l', 's', 'Fo', 'Fs', 'dummy']
    elif ndat == 6:
        names = ['h', 'k', 'l', 's', 'Fo', 'Fs']
    else:
        raise ValueError(f'Unexpected number of data columns: {ndat}')
    hkl = pd.read_csv(fpath, sep=r'\s+', skiprows=1, names=names)
    return head, hkl, ndat

def cart_from_cell(cell: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Converts crystallographic unit cell parameters to Cartesian lattice vectors.

    Parameters
    ----------
    cell : np.ndarray
        A 1D numpy array of shape (6,) containing the unit cell parameters:
        [a, b, c, alpha, beta, gamma], where a, b, c are the cell lengths in Ångströms,
        and alpha, beta, gamma are the cell angles in degrees.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple containing three numpy arrays representing the Cartesian lattice vectors
        (a, b, c), each of shape (3,).

    Raises
    ------
    ValueError
        If the input array does not have shape (6,).

    Notes
    -----
    The returned vectors are in meters (converted from Ångströms).
    """
    if cell.shape != (6,):
        raise ValueError('Lattice constants must be 1d array with 6 elements')
    a, b, c = cell[:3]*1E-10
    alpha, beta, gamma = np.radians(cell[3:])
    av = np.array([a, 0, 0], dtype=float)
    bv = np.array([b * np.cos(gamma), b * np.sin(gamma), 0], dtype=float)
    # calculate vector c
    x = np.cos(beta)
    y = (np.cos(alpha) - x * np.cos(gamma)) / np.sin(gamma)
    z = np.sqrt(1. - x**2. - y**2.)
    cv = np.array([x, y, z], dtype=float)
    cv /= np.linalg.norm(cv)
    cv *= c
    return av, bv, cv

def matrix_from_cell(cell):
    """
    Computes the reciprocal lattice transformation matrix from unit cell parameters.

    Given a unit cell described by its parameters, this function calculates the transformation
    matrix (A) that converts real-space vectors to reciprocal-space vectors. The matrix is
    constructed using the reciprocal lattice vectors a*, b*, and c*, which are derived from
    the real-space vectors a, b, and c.

    Parameters:
        cell (array-like): A sequence or array containing the unit cell parameters,
            typically [a, b, c, alpha, beta, gamma] or equivalent.

    Returns:
        numpy.ndarray: A 3x3 matrix representing the transformation from real-space
            to reciprocal-space, rounded to 6 decimal places.

    Notes:
        - The function assumes the existence of a helper function `cart_from_cell(cell)`
          that returns the real-space vectors (a, b, c) in Cartesian coordinates.
        - The returned matrix columns correspond to the reciprocal lattice vectors a*, b*, and c*.
    """
    cell = np.array(cell)
    av, bv, cv = cart_from_cell(cell)
    a_star = (np.cross(bv, cv)) / (np.cross(bv, cv).dot(av))
    b_star = (np.cross(cv, av)) / (np.cross(cv, av).dot(bv))
    c_star = (np.cross(av, bv)) / (np.cross(av, bv).dot(cv))
    A = np.zeros((3, 3), dtype='float')  # transform matrix
    A[:, 0] = a_star
    A[:, 1] = b_star
    A[:, 2] = c_star
    return np.round(A, 6) * 1e-10

def write_hkl(head, data, ndat, fpath):
    """
    Writes reflection data to an HKL file with a specified header.

    Parameters:
        head (str): The header string to write at the top of the file.
        data (pandas.DataFrame): DataFrame containing columns 'h', 'k', 'l', 's', 'Fo', and 'Fs'.
        fpath (str): The file path where the HKL file will be written.

    The function writes the header followed by each row of the DataFrame in a formatted manner.
    Each row is written with fixed-width columns for h, k, l, s, Fo, and Fs.
    """
    with open(fpath, 'w') as f:
        f.write(head)
        #   0   0   4 1  142.372624    0.337509
        if ndat == 6:
            for row in data.itertuples(index=False):
                f.write(f'{row.h:4}{row.k:4}{row.l:4}{row.s:4}{row.Fo:12}{row.Fs:12}\n')
        elif ndat == 7:
            for row in data.itertuples(index=False):
                f.write(f'{row.h:4}{row.k:4}{row.l:4}{row.s:4}{row.Fo:12}{row.Fs:12}{row.dummy:8}\n')
        else:
            raise ValueError(f'Unexpected number of data columns: {ndat}')
    print(f'Data written to {fpath}')

def main():
    # Check if the current working directory contains 'xd.mas'
    homedir = os.path.dirname(__file__)
    if os.path.exists(os.path.join(homedir, 'xd.mas')): pass
    else: raise SystemExit(f'xd.mas not found in {homedir}')
    # Define the number of scale factors
    sf_num = 10
    # Get the cell parameters and wavelength from the XD file
    cell, wave = get_cell_and_wave(os.path.join(homedir, 'xd.mas'))
    # Get the header and data from the HKL file
    # Data format: DataFrame containing columns ['h', 'k', 'l', 's', 'Fo', 'Fs']
    head, data, ndat = get_data(os.path.join(homedir, 'xd.hkl'))
    # Convert the cell parameters to a transformation matrix
    matrix = matrix_from_cell(cell)
    # Prepare the data for scaling
    hkl = data[['h', 'k', 'l']].values
    # Calculate the scattering vector in reciprocal space
    scatvec = matrix.dot(hkl.T).T
    # Calculate the resolution of the scattering vector as sin(theta)/lambda
    # The resolution is calculated as the norm of the scattering vector.
    stl = np.linalg.norm(scatvec, axis=1) / 2
    # Create 'scalefactor' number of resolution bins
    sf_range = np.linspace(0, stl.max(), sf_num + 1)
    # Check if the xd.fco file exists and if its stl values 
    # are consistent with the calculated stl values
    if os.path.exists(os.path.join(homedir, 'xd.fco')):
        fco = pd.read_csv(os.path.join(homedir, 'xd.fco'), sep=r'\s+', skiprows=26, names=['h', 'k', 'l', 'Fc', 'Fo', 'Fs', 'stl', 'dummy'])
        assert np.isclose(fco['stl'].max(), sf_range.max(), atol=1E-4), 'Calculation inconsistent with the xd.fco stl values!'
    else:
        print('xd.fco not found, skipping consistency check!')
    # Assign scale factors based on the resolution bins
    data['s'] = pd.cut(stl, bins=sf_range, labels=np.arange(1, sf_num+1), include_lowest=True)
    # Write the scaled data to a new HKL file
    write_hkl(head, data, ndat, os.path.join(homedir, 'xd_scaled.hkl'))

if __name__ == '__main__':
    main()