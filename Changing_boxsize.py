from simassis.files import read_configuration, configuration
import numpy as np
import shutil
import os


def get_basis(simassis_configuration):
    basis_M = np.array([simassis_configuration.basex,
                        simassis_configuration.basey,
                        simassis_configuration.basez])
    return basis_M

def get_cartesian_coordinates(cell_splited, new_base, atoms_variety):
    cartesian_coordinates_splited = []
    string = ""
    for particle in range(atoms_variety):
        cartesian_coordinates_splited.append(
            np.dot(new_base, cell_splited[particle].T).T)
    return cartesian_coordinates_splited

def add_unit_cell_to_shape(shape_array, simassis_configuration):
    x = shape_array[0]
    y = shape_array[1]
    z = shape_array[2]
    basis_M = get_basis(simassis_configuration)
    extended_cell_splited = []
    for particle in range(simassis_configuration.atoms_variety):
        extended_cell_splited.append([])

    new_base_vec = np.array([x, y, z])
    new_base_vec = np.dot(basis_M, new_base_vec)

    for dx in range(x):
        for dy in range(y):
            for dz in range(z):
                move_vec = np.array([dx, dy, dz])
                move_vec = np.dot(basis_M, move_vec.T).T
                #print(move_vec)
                for particle in range(simassis_configuration.atoms_variety):
                    normalaised_result = \
                        (simassis_configuration.positions_splited[particle] 
                        + move_vec) \
                        / new_base_vec
                    extended_cell_splited[particle].append(normalaised_result)

    for particle in range(simassis_configuration.atoms_variety):           
        extended_cell_splited[particle] = \
            np.concatenate(extended_cell_splited[particle])

    new_base = np.array([[new_base_vec[0], 0, 0],
        [0, new_base_vec[1], 0],
        [0, 0, new_base_vec[2]]])

    return new_base, extended_cell_splited

def new_atoms_quantity(new_dimensions, simassis_configuration):
    new_atoms_quantity = []
    multiplier = new_dimensions.prod()

    for particle in range(simassis_configuration.atoms_variety):  
        new_atoms_quantity.append(
            simassis_configuration.atoms_quantity[particle] * multiplier)

    return new_atoms_quantity

def create_extended_cell(new_dimensions, 
    path = "POSCAR.vasp", 
    output_filename = "POSCAR.vasp"):
    conf = read_configuration(path, "vasp")
    new_dimensions = np.array(new_dimensions)
    new_base, extended_cell_splited = add_unit_cell_to_shape(new_dimensions, conf)
    new_conf = configuration(
        name = "extended_box",
        basex = new_base[0,:],
        basey = new_base[1,:],
        basez = new_base[2,:],
        cart_dir = "Direct",
        atoms_type = conf.atoms_type,
        atoms_quantity = new_atoms_quantity(new_dimensions, conf),
        atoms_variety = conf.atoms_variety,
        positions_splited = extended_cell_splited
        )
    new_conf.write('vasp', output_filename)
    print("New File: " + output_filename)

def create_vacancies( 
    vacancy_procentage = 3, 
    input_filename = "POSCAR.vasp",
    output_filename = "POSCAR.vasp"):

    f = open(input_filename, "r")
    f_lines = f.readlines()
    f_length = len(f_lines)
    f.close()
    amount_of_oxygenes = int(f_lines[6].split()[1])
    print(amount_of_oxygenes)
    amount_of_vacanies = int(amount_of_oxygenes/100 * vacancy_procentage)

    rng = np.random.default_rng()
    random_indexes = rng.choice(
        amount_of_oxygenes - 1 - amount_of_vacanies, 
        size = amount_of_vacanies , 
        replace = False) + 1 # Replace flag makes all numbers unique, +1 we dont want 0
    print(random_indexes)

    for i in random_indexes: 
        f_lines.pop(-i)

    f_lines[6] = '\t' + f_lines[6].split()[0] \
        + '\t' + str( amount_of_oxygenes - amount_of_vacanies ) + '\n'

    f = open("POSCAR.vasp", "w")
    f.writelines(f_lines)
    f.close()

    print("New POSCAR was created without " 
        + str(amount_of_vacanies) 
        + "-O (" + str(vacancy_procentage) +"%)")
    #print("Deleted lines: " + str(f_length - random_indexes))

def convert_to_lammps_data(input_filename = "POSCAR.vasp"):
    conf = read_configuration(input_filename, "vasp")
    conf.write('lammps', 'data', ['triclinic', 'charge_custom',[0,0]])

def create_local_copy_of_POSCAR(path):
    shutil.copy(path, ".")

if __name__ == '__main__':
    path_to_POSCAR = r"/home/tymon/Desktop/Lammps/Zmiana rozmiaru CeO/src/POSCAR.vasp"

    create_local_copy_of_POSCAR(path_to_POSCAR)
    create_extended_cell([25, 25, 25])
    create_vacancies(3)
    convert_to_lammps_data()