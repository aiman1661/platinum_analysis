import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pandas import DataFrame as df
import dataframe_image as dfi

pd.set_option('expand_frame_repr', False)

sin = np.sin
rad = np.radians
sqrt = np.sqrt

def sqrt_ind(miller_list):
    "Calculate the sqrt of square sum of miller indices and put in dictionary"
    sqrt_sum = []
    for hkl in miller_list:
        g_hkl = sqrt(hkl[0]**2 + hkl[1]**2 + hkl[2]**2)
        sqrt_sum.append(g_hkl)

    return sqrt_sum

def hkl_generator(index):
    "Calculate possible miller indices (hkl and order"
    miller_list = []
    for i in index:
        h = i
        for j in index:
            k = j 
            for m in index:
                l = m
                miller_ind = [h,k,l]
                miller_list.append(miller_ind)

    for plane in miller_list:
        plane.sort()
    
    res_miller_list = [] # Unique (hkl) up to permutation and negatives

    for plane in miller_list:
        if plane not in res_miller_list:
            res_miller_list.append(plane)

    return res_miller_list[1:]

def fcc_indices(miller_list):
    "Listing possible fcc miller indices"
    final_list = []

    for plane in miller_list:
        odd_counter = 0
        even_counter = 0
        for i in plane:
            if (i % 2) == 0:
                even_counter += 1
            else:
                odd_counter += 1
        
        if even_counter == 3 or odd_counter == 3:
            final_list.append(plane)

    return final_list

def bcc_indices(miller_list):
    "Listing possible bcc miller indices"
    final_list = []

    for plane in miller_list:
        if sum(plane) % 2 == 0:
            final_list.append(plane)

    return final_list

def lattice_param_a(d_space, hkl_list):
    "Calculate values of lattice parameter using hkl"
    param_a_list = []
    for plane in hkl_list:
        a = d_space * sqrt(plane[0]**2 + plane[1]**2 + plane[2]**2)
        param_a_list.append(np.array(a))

    return np.array(param_a_list)

def plotter(d_space, miller_ind, a_params, title):
    "Plot the data for either bcc or fcc, and create dataframes"

    y_axis_data = a_params

    fig, ax = plt.subplots(figsize=(10,5))

    x_axis_labels = [f'({x[0]}, {x[1]}, {x[2]})' for x in miller_ind]

    dataframe_list = []
    g_hkl = sqrt_ind(miller_ind)
    dataframe_g_hkl = df(data = g_hkl, index = x_axis_labels, columns = ["g_hkl"])

    for i in range(a_params.shape[1]):
        y_axis = y_axis_data[:,i]
        plt.scatter(x_axis_labels,y_axis, label=f'd={d_space[i]}', marker='x')

        dataframe = df(data = y_axis, index = x_axis_labels, columns = [f'd={d_space[i]}'])
        dataframe_list.append(dataframe)

    ax.set_xlabel('(h k l)')
    ax.set_ylabel(r'Lattice parameter, a ($\AA$)')
    ax.set_title(f'{title}')
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5), fancybox=True, shadow=True)
    plt.xticks(rotation=45)  # Rotate x-axis labels for clarity
    plt.tight_layout()
    plt.show()

    full_df_carry = dataframe_g_hkl.join(dataframe_list[:])
    sorted_df_carry = full_df_carry.sort_values(by = 'g_hkl')
    sorted_df_red_carry = sorted_df_carry.drop(sorted_df_carry.index[9:]) # Reduced (to first 9) sorted df

    dfi.export(sorted_df_red_carry,f'{title}.png')
    sorted_df_red_carry.to_excel(f"{title}_out.xlsx")
    return print(sorted_df_red_carry)

def main():

    # Data and parameter given
    diff_angle = np.array([10.378, 11.992, 16.999, 19.955, 20.858, 24.118, 26.332, 27.019, 29.661]) # In degrees, 2theta
    theta = diff_angle/2.0
    wavelength = 0.41000 # In Angstrom (E-10)

    # Range of possible miller integers
    index_set = np.arange(5)

    # Use Bragg's law to calculate d-spacing
    d_space = wavelength/(2*sin(rad(theta))) # In Angstrom

    d_space_df = df({"Diffraction angle": diff_angle, 
                     "d_hkl": d_space})
    print(d_space_df)
    dfi.export(d_space_df,'d_space.png')
    print('\n')

    # List of all possible miller indices
    miller_list = hkl_generator(index_set)

    # Lists of possible miller indices for fcc and bcc, taking systematic absences into account
    miller_fcc = fcc_indices(miller_list)
    miller_bcc = bcc_indices(miller_list)

    # Lattice parameter 'a' values
    lattice_param_fcc = lattice_param_a(d_space, miller_fcc)
    lattice_param_bcc = lattice_param_a(d_space, miller_bcc)

    print('FCC Measurements')
    plotter(d_space,miller_fcc,lattice_param_fcc, 'FCC Measurements')
    print('\n')
    print('BCC Measurements')
    plotter(d_space,miller_bcc,lattice_param_bcc, 'BCC Measurements')
    print('\n')
    

if __name__ == '__main__':
    main()
    