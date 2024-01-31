import numpy as np
import sys
 
def rotation_matrix(axis, angle):

    # Convert the axis to a unit vector
    axis = axis / np.linalg.norm(axis)

    # Convert the angle from degrees to radians
    angle_rad = np.deg2rad(angle)

    # Compute the cross product matrix of the rotation axis
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])

    # Compute the rotation matrix using Rodrigues' rotation formula
    R = np.eye(3) + np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * np.dot(K, K)
    R = np.round(R, 3)
    return R.T

def mirror_image(axis):

    # Convert the reflection axis to a unit vector
    axis = axis / np.linalg.norm(axis)

    
    # Compute the reflection matrix
    R_reflection = np.eye(3) - 2 * np.outer(axis, axis)


    # Combine the rotation and reflection matrices
    R_mirror = R_reflection
    return R_mirror

def Improper_rotation(axis,angle):
    
    # Convert the axis to a unit vector
    axis = axis / np.linalg.norm(axis)

    # Convert the angle from degrees to radians
    angle_rad = np.deg2rad(angle)

    # Compute the cross product matrix of the rotation axis
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])

    # Compute the rotation matrix using Rodrigues' rotation formula
    R = np.eye(3) + np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * np.dot(K, K)
    R = np.round(R, 3)
    
    # Convert the reflection axis to a unit vector
    axis = axis / np.linalg.norm(axis)

    
    # Compute the reflection matrix
    R_reflection = np.eye(3) - 2 * np.outer(axis, axis)


    # Combine the rotation and reflection matrices
    R_Improper = np.dot(R.T,R_reflection)
    return R_Improper


def test_pointgroup(Operations,Structure,Inp_Atoms,dev_name):
    syms = []
    for index,operation in enumerate(Operations):
        operated_structure =  np.dot(Structure, operation)
        operated_structure_permuted,operated_atoms_permuted,permutation = best_permutation_multiple_atoms(Structure,Inp_Atoms,operated_structure)
        s = symmetry_deviation_calculator(operated_structure_permuted, Structure)
        syms.append(s)
        print('deviation',dev_name[index],s)
    print(np.mean(np.array(syms)))
    return np.array(syms)

def SaveSymmetryFile(PointGroup,OperationNames,Operations):
    # Saving the lists to a CSV file
    data = list(zip(Operations,OperationNames))
    filename = 'pointgroup_{}.csv'.format(PointGroup)
    
    dirname = os.path.dirname(__file__)
    Path_To_Symmetry_Operations = os.path.join(dirname, 'Operations')
    
    save_folderpath = Path_To_Symmetry_Operations
    if not os.path.exists(save_folderpath):
        os.makedirs(save_folderpath)
                    
    with open(os.path.join(save_folderpath,filename), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['OperationName', 'Operation'])  # Write header
        for matrix,name in data:
            matrix = np.round(matrix,4)
            matrix = matrix.reshape(1,9)[0]
            writer.writerow([name,matrix[0],matrix[1],matrix[2],matrix[3],matrix[4],matrix[5],matrix[6],matrix[7],matrix[8]])
    
    print(f'Saved data to {filename}')

def ReadSymmetryFile(PointGroup):
    # Reading the lists from the CSV file
    read_name_list = []
    read_matrix_list = []
    
    dirname = os.path.dirname(__file__)
    Path_to_folder = os.path.join(dirname, 'Operations')
    
    with open(os.path.join(Path_to_folder,PointGroup + '.csv'), 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            name = row[0]
            matrix_str = [row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9]]
            matrix_str = np.array(matrix_str, dtype=float).reshape(3,3)
    
            read_name_list.append(name)
            read_matrix_list.append(matrix_str)
    
    print('Read data from CSV:', PointGroup)
    return read_name_list,read_matrix_list

#%%
# =============================================================================
# =============================================================================
# # Operations
# =============================================================================
# =============================================================================

Symmetry_Operations = []
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]))) #i
Save_And_Read_SymOperations.SaveSymmetryFile('operation_i',point_group_symmetry_names,Symmetry_Operations)


Symmetry_Operations = []
axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # sigmav
point_group_symmetry_names = ['sigma_v']
Save_And_Read_SymOperations.SaveSymmetryFile('operation_sigma_v',point_group_symmetry_names,Symmetry_Operations)


Symmetry_Operations = []
axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigmav
point_group_symmetry_names = ['sigma_v']
Save_And_Read_SymOperations.SaveSymmetryFile('operation_sigma_h',point_group_symmetry_names,Symmetry_Operations)



Symmetry_Operations = []
axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
point_group_symmetry_names = ['C2']
Save_And_Read_SymOperations.SaveSymmetryFile('operation_C2',point_group_symmetry_names,Symmetry_Operations)

Symmetry_Operations = []
axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C3
point_group_symmetry_names = ['C3']
Save_And_Read_SymOperations.SaveSymmetryFile('operation_C3',point_group_symmetry_names,Symmetry_Operations)

Symmetry_Operations = []
axis = [0,0,1]
angle = 240
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C3_2
point_group_symmetry_names = ['C3_2']
Save_And_Read_SymOperations.SaveSymmetryFile('operation_C3_2',point_group_symmetry_names,Symmetry_Operations)

Symmetry_Operations = []
axis = [0,0,1]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # S4
point_group_symmetry_names = ['S4']
Save_And_Read_SymOperations.SaveSymmetryFile('operation_S4',point_group_symmetry_names,Symmetry_Operations)





#%%
# =============================================================================
# =============================================================================
# # pointgroup_C2
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 90*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2

point_group_symmetry_names = ['C2']
Save_And_Read_SymOperations.SaveSymmetryFile('C2',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C2')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_C3
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C3
axis = [0,0,1]
angle = 120*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C3^2


point_group_symmetry_names = ['C3','C3^2']
Save_And_Read_SymOperations.SaveSymmetryFile('C3',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C3')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_C4
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C4
axis = [0,0,1]
angle = 90*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
axis = [0,0,1]
angle = 90 * 3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C4^3


point_group_symmetry_names = ['C4', 'C2', 'C4^3']
Save_And_Read_SymOperations.SaveSymmetryFile('C4',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C4')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_D4d
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 45
Symmetry_Operations.append(Improper_rotation(axis,angle)) # 2S8
axis = [0,0,1]
angle = 45*7
Symmetry_Operations.append(Improper_rotation(axis,angle)) # 2S8
axis = [0,0,1]
angle = 90

Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
axis = [0,0,1]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
axis = [0,0,1]
angle = 45 *3

Symmetry_Operations.append(Improper_rotation(axis,angle)) # 2(S8)3
axis = [0,0,1]
angle = 45*5
Symmetry_Operations.append(Improper_rotation(axis,angle)) # 2(S8)3
axis = [0,0,1]
angle = 90*2

Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
axis = [np.cos(np.pi/180 * 22.5),np.cos(np.pi/180 * 67.5),0]
angle = 180

Symmetry_Operations.append(rotation_matrix(axis, angle)) # 4C'2
axis = [np.cos(np.pi/180 *( 22.5+90)),np.cos(np.pi/180 * (67.5 +90)),0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 4C'2
axis = [np.cos(np.pi/180 *( 22.5+180)),np.cos(np.pi/180 * (67.5 +180)),0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 4C'2
axis = [np.cos(np.pi/180 *( 22.5+3*90)),np.cos(np.pi/180 * (67.5 +3*90)),0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 4C'2

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # 4sigma_d
axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 4sigma_d
axis =[1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 4sigma_d
axis =[-1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 4sigma_d

point_group_symmetry_names = ['2S8','2S8','2C4','2C4','2(S8)3','2(S8)3','C2','4C^2','4C^2','4C^2','4C^2','4sigma_d','4sigma_d','4sigma_d','4sigma_d']
Save_And_Read_SymOperations.SaveSymmetryFile('D4d',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D4d')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_D3h
# =============================================================================
# =============================================================================

Symmetry_Operations = []

axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3
angle = 120 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C3


axis = [-np.sqrt(3)/2,1/2,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis =[np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 3C2
axis =[0,-1,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h


axis = [0,0,1]
angle = 120
Symmetry_Operations.append( Improper_rotation(axis,angle)) # 2S3
angle = 120 *2 
Symmetry_Operations.append( Improper_rotation(axis,angle)) # 2S3



axis =[0,0,1]
Symmetry_Operations.append( mirror_image(axis)) # 3sigma_v
axis =[1/2,np.sqrt(3)/2,0]
Symmetry_Operations.append( mirror_image(axis)) # 3sigma_v
axis = [1/2,-np.sqrt(3)/2,0]
Symmetry_Operations.append( mirror_image(axis)) # 3sigma_v

point_group_symmetry_names = ['2C3','2C3','3C2','3C2','3C2','sigma_h','2S3','2S3','3sigma_v','3sigma_v','3sigma_v']
Save_And_Read_SymOperations.SaveSymmetryFile('D3h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D3h')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_C4v
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
axis = [0,0,1]
angle = 90 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C4

axis = [0,0,1]
angle = 90 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 1C2

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # 2sigma_v

axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 2sigma_v

axis =[1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 2sigma_d

axis =[-1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 2sigma_d

point_group_symmetry_names = ['2C4','2C4','1C2','2sigma_v','2sigma_v','2sigma_d','2sigma_d']
Save_And_Read_SymOperations.SaveSymmetryFile('C4v',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C4v')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_C5v
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 360 / 5
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2
axis = [0,0,1]
angle = 360 / 5 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis =[np.cos(2*np.pi / 5),np.sin(2*np.pi / 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis =[np.cos(2*np.pi / 5*2),np.sin(2*np.pi / 5 *2),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis =[np.cos(2*np.pi / 5*3),np.sin(2*np.pi / 5 * 3),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis =[np.cos(2*np.pi / 5*4),np.sin(2*np.pi / 5*4),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v



point_group_symmetry_names = ['2C5','2C5',
                              '2C5^2','2C5^2',
                              '5sigma_v','5sigma_v','5sigma_v','5sigma_v','5sigma_v'
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('C5v',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C5v')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_C5
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 360 / 5
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2
axis = [0,0,1]
angle = 360 / 5 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2


point_group_symmetry_names = ['2C5','2C5',
                              '2C5^2','2C5^2',
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('C5',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C5')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_C3v
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3
axis = [0,0,1]
angle = 120 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C3

axis = [-np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v

axis =[np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append( mirror_image(axis)) # 3sigma_v
axis =[0,-1,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v


point_group_symmetry_names = ['2C3','2C3'
                              '3sigma_v','3sigma_v','3sigma_v']

Save_And_Read_SymOperations.SaveSymmetryFile('C3v',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C3v')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D6h
# =============================================================================
# =============================================================================
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
angle = 360 / 2 *5
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle =  360 / 2 *2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6^2
angle = 360 / 2 *4
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6^2
# =============================================================================
# C2
# =============================================================================
angle = 360 / 2 *3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
# =============================================================================
# 3C_2
# =============================================================================
angle = 180
axis = [1,0,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C_2
axis = [np.cos(np.pi*2 / 6),np.sin(np.pi*2 / 6),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C_2
axis = [np.cos(np.pi*2 / 6 * 2),np.sin(np.pi*2 / 6 * 2),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C_2


# =============================================================================
# 3C__2
# =============================================================================
angle = 180
axis = [np.cos(np.pi*2 / 12),np.sin(np.pi*2 / 12),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C__2
axis = [np.cos(np.pi*2 / 12 * 3),np.sin(np.pi*2 / 12* 3),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C__2
axis = [np.cos(np.pi*2 / 12 * 5),np.sin(np.pi*2 / 12* 5),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C__2
# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]))) #i

# =============================================================================
# 2S3
# =============================================================================
axis = [0,0,1]
angle = 360 / 3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S3
angle = 360 / 3 *2
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S3
# =============================================================================
# 2S6
# =============================================================================
axis = [0,0,1]
angle =  360 / 6
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S6
angle = 360 / 6 *5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S6



# =============================================================================
#  3sigma_h
# =============================================================================
axis = [0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_h
# =============================================================================
#  3sigma_v
# =============================================================================
axis = [1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v
axis = [np.cos(np.pi*2 / 6),np.sin(np.pi*2 / 6),0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v
axis = [np.cos(np.pi*2 / 6 * 2),np.sin(np.pi*2 / 6 * 2),0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v
# =============================================================================
#  3sigma_d
# =============================================================================
axis = [np.cos(np.pi*2 / 12),np.sin(np.pi*2 / 12),0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d
axis = [np.cos(np.pi*2 / 12 * 3),np.sin(np.pi*2 / 12* 3),0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d
axis = [np.cos(np.pi*2 / 12 * 5),np.sin(np.pi*2 / 12* 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d
point_group_symmetry_names = ['2C6','2C6',
                              '2C6^2','2C6^2',
                              'C2',
                              '3C_2','3C_2','3C_2',
                              '3C__2','3C__2','3C__2',
                              'i',
                              '2S3','2S3',
                              '2S6','2S6',
                              '3sigma_h',
                              '3sigma_v','3sigma_v','3sigma_v',
                              '3sigma_d','3sigma_d','3sigma_d',

                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('D6h',point_group_symmetry_names,Symmetry_Operations)
symname,sumop = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D6h')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_S6
# =============================================================================
# =============================================================================
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
angle = 360 / 3 *2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6

# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]))) #i

# =============================================================================
# 2S6
# =============================================================================
axis = [0,0,1]
angle =  360 / 6
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S6
angle = 360 / 6 *5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S6

point_group_symmetry_names = ['2C6','2C6',
                              'i',
                              '2S6','2S6',
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('S6',point_group_symmetry_names,Symmetry_Operations)
symname,sumop = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_S6')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D6
# =============================================================================
# =============================================================================
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
angle = 360 / 2 *5
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle =  360 / 2 *2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6^2
angle = 360 / 2 *4
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6^2
# =============================================================================
# C2
# =============================================================================
angle = 360 / 2 *3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
# =============================================================================
# 3C_2
# =============================================================================
angle = 180
axis = [1,0,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C_2
axis = [np.cos(np.pi*2 / 6),np.sin(np.pi*2 / 6),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C_2
axis = [np.cos(np.pi*2 / 6 * 2),np.sin(np.pi*2 / 6 * 2),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C_2


# =============================================================================
# 3C__2
# =============================================================================
angle = 180
axis = [np.cos(np.pi*2 / 12),np.sin(np.pi*2 / 12),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C__2
axis = [np.cos(np.pi*2 / 12 * 3),np.sin(np.pi*2 / 12* 3),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C__2
axis = [np.cos(np.pi*2 / 12 * 5),np.sin(np.pi*2 / 12* 5),0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C__2

point_group_symmetry_names = ['2C6','2C6',
                              '2C6^2','2C6^2',
                              'C2',
                              '3C_2','3C_2','3C_2',
                              '3C__2','3C__2','3C__2',
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('D6',point_group_symmetry_names,Symmetry_Operations)
symname,sumop = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D6')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_C6
# =============================================================================
# =============================================================================
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 6
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
angle = 360 / 6 *5
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6
# =============================================================================
# 2C6
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle =  360 / 6 *2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6^2
angle = 360 / 6 *4
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C6^2
# =============================================================================
# C2
# =============================================================================
angle = 360 / 6 *3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2

point_group_symmetry_names = ['2C6','2C6',
                              '2C6^2','2C6^2',
                              'C2',
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('C6',point_group_symmetry_names,Symmetry_Operations)
symname,sumop = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C6')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_D3d
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3
angle = 120*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3


angle = 180
axis = [-np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis =[np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 3C2
axis =[0,-1,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]))) #i
# =============================================================================
# 2S6
# =============================================================================
axis = [0,0,1]
angle = 120/2
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2C3
angle = 120/2 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2C3
# =============================================================================
#  3sigma_d
# =============================================================================
axis = [-np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d
axis =[np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append( mirror_image(axis)) # 3sigma_d
axis =[0,-1,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d


point_group_symmetry_names = ['2C3','2C3',
                              '3C2','3C2','3C2',
                              'i',
                              '2S6','2S6',
                              '3sigma_d','3sigma_d','3sigma_d']

Save_And_Read_SymOperations.SaveSymmetryFile('D3d',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D3d')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D4h
# =============================================================================
# =============================================================================
# =============================================================================
# 2C4
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
# =============================================================================
# C2
# =============================================================================
angle = 90*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
# =============================================================================
# 2C_2
# =============================================================================
angle = 180
axis = [1,0,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C_2
axis = [0,1,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C_2
# =============================================================================
# 2C__2
# =============================================================================
axis = [1,1,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C__2
axis = [1,-1,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C__2
# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]))) #i
# =============================================================================
# 2S4
# =============================================================================
axis = [0,0,1]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S4
angle = 90*3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S4

# =============================================================================
#  3sigma_h
# =============================================================================
axis = [0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_h
# =============================================================================
#  3sigma_v
# =============================================================================
axis = [1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v
axis = [0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_v
# =============================================================================
#  3sigma_d
# =============================================================================
axis = [1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d
axis = [1,-1,0]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_d


point_group_symmetry_names = ['2C4','2C4',
                              'C2',
                              'C2_',
                              'C2__','C2__',
                              'i',
                              '2S4','2S4',
                              'sigma_h',
                              '2sigma_v','sigma_v',
                              '2sigma_d','sigma_d'
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('D4h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D4h')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_C4h
# =============================================================================
# =============================================================================
# =============================================================================
# 2C4
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
# =============================================================================
# C2
# =============================================================================
angle = 90*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2

# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]))) #i
# =============================================================================
# 2S4
# =============================================================================
axis = [0,0,1]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S4
angle = 90*3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S4

# =============================================================================
#  3sigma_h
# =============================================================================
axis = [0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_h


point_group_symmetry_names = ['2C4','2C4',
                              'C2',
                              'i',
                              '2S4','2S4',
                              'sigma_h',
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('C4h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C4h')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_C3h
# =============================================================================
# =============================================================================
# =============================================================================
# 2C3
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3
angle = 120*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3


# =============================================================================
# 2S3
# =============================================================================
axis = [0,0,1]
angle = 120
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S3
angle = 120*2
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S3

# =============================================================================
#  3sigma_h
# =============================================================================
axis = [0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # 3sigma_h


point_group_symmetry_names = ['2C3','2C3',
                              '2S3','2S3',
                              'sigma_h',
                              ]

Save_And_Read_SymOperations.SaveSymmetryFile('C3h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C3h')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_Cs
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigmah


point_group_symmetry_names = ['sigma_h']

Save_And_Read_SymOperations.SaveSymmetryFile('Cs',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_Cs')

#%%
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D2d
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 2S4
angle = 90 *3
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S4

axis = [0,0,1]
angle = 90 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 1C2

axis = [1,1,0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C2
axis = [1,-1,0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C2

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # 2sigma_d
axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 2sigma_d

point_group_symmetry_names = ['2S4','2S4','1C2', '2C2', '2C2', '2sigma_d','2sigma_d']
Save_And_Read_SymOperations.SaveSymmetryFile('D2d',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D2d')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_Oh
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 8C3
# =============================================================================
axis = [1,1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

# =============================================================================
# 6C2
# =============================================================================
axis = [1,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [1,0,-1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [0,1,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [0,1,-1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [1,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [1,-1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

# =============================================================================
# 6C4
# =============================================================================
axis = [1,0,0]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4
axis = [1,0,0]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4

axis = [0,1,0]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4
axis = [0,1,0]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4
axis = [0,0,1]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4

# =============================================================================
# 3C2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]])))

# =============================================================================
# 6S4
# =============================================================================

axis = [1,0,0]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4
axis = [1,0,0]
angle = 90*3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

axis = [0,1,0]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4
axis = [0,1,0]
angle = 90*3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4
axis = [0,0,1]
angle = 90*3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

# =============================================================================
# 8S6
# =============================================================================
axis = [1,1,1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6
axis = [1,1,1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,1,-1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,1,-1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,-1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,-1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

# =============================================================================
# 3sigma_h
# =============================================================================

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h

axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h
axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h

# =============================================================================
# 6sigma_d
# =============================================================================
axis =[1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[1,0,1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[0,1,1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[1,-1,0]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[1,0,-1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[0,1,-1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d

point_group_symmetry_names = ['8C3','8C3','8C3','8C3','8C3','8C3','8C3','8C3',
                              '6C2','6C2','6C2','6C2','6C2','6C2',
                              '6C4','6C4','6C4','6C4','6C4','6C4',
                              '3C2','3C2','3C2',
                              'i',
                              '6S4','6S4','6S4','6S4','6S4','6S4',
                              '8S6','8S6','8S6','8S6','8S6','8S6',
                              'sigma_h','sigma_h','sigma_h',
                              '6sigma_d','6sigma_d','6sigma_d','6sigma_d','6sigma_d','6sigma_d','6sigma_d','6sigma_d'                              
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('Oh',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_Oh')[1]


#%%
# =============================================================================
# =============================================================================
# # pointgroup_O
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 8C3
# =============================================================================
axis = [1,1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

# =============================================================================
# 6C2
# =============================================================================
axis = [1,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [1,0,-1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [0,1,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [0,1,-1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [1,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

axis = [1,-1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C2

# =============================================================================
# 6C4
# =============================================================================
axis = [1,0,0]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4
axis = [1,0,0]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4

axis = [0,1,0]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4
axis = [0,1,0]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4
axis = [0,0,1]
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 6C4

# =============================================================================
# 3C2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

point_group_symmetry_names = ['8C3','8C3','8C3','8C3','8C3','8C3','8C3','8C3',
                              '6C2','6C2','6C2','6C2','6C2','6C2',
                              '6C4','6C4','6C4','6C4','6C4','6C4',
                              '3C2','3C2','3C2',                        
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('O',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_O')[1]

Inp_Atoms,Structure = load_xyz(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\ModelPolyhedron_8-9-10\CN8_Square.xyz',2)
test_pointgroup(Operations,Structure,Inp_Atoms)
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D4
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 2C4
# =============================================================================
axis = [0,0,1]
angle = 90
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
angle = 90*3
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C4
# =============================================================================
# C2
# =============================================================================
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
# =============================================================================
# 2C'2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C'2
axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C'2
# =============================================================================
# 2C''2
# =============================================================================
axis = [1,-1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C"2
axis = [1,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C"2

point_group_symmetry_names = ['2C4','2C4',
                              '2C2',
                              '2C2','2C2',
                              '2C2','2C2'
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('D4',point_group_symmetry_names,Symmetry_Operations)
dev_name, Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D4')

Inp_Atoms,Structure = load_xyz(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\Structures_in_different_point_groups\pointgroup_D4d_S8.xyz',1)

test_pointgroup(Operations,Structure,Inp_Atoms,dev_name)
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D3
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 2C3
# =============================================================================
axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3
angle = 120*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 2C3
# =============================================================================
# C2
# =============================================================================
axis = [-np.sqrt(3)/2,1/2,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis =[np.sqrt(3)/2,1/2,0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 3C2
axis =[0,-1,0]
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2



point_group_symmetry_names = ['2C3','2C3',
                              '3C2','3C2','3C2'

                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('D3',point_group_symmetry_names,Symmetry_Operations)
dev_name, Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D3')

Inp_Atoms,Structure = load_xyz(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\Structures_in_different_point_groups\pointgroup_D3_[Co(C2H8N2)3]Cl2.xyz',1)

test_pointgroup(Operations,Structure,Inp_Atoms,dev_name)


#%%
# =============================================================================
# =============================================================================
# # pointgroup_Td
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 8C3
# =============================================================================
axis = [1,1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

# =============================================================================
# 3C2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

# =============================================================================
# 6S4
# =============================================================================
axis = [1,0,0]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

axis = [0,1,0]
angle = 90 
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

axis = [0,0,1]
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

axis = [1,0,0]
angle = -90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4


axis = [0,1,0]
angle = -90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

axis = [0,0,1]
angle = -90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 6S4

# =============================================================================
# 6sigma_d
# =============================================================================
axis =[1,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[1,0,1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[0,1,1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[1,-1,0]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[1,0,-1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d
axis =[0,1,-1]
Symmetry_Operations.append(mirror_image(axis)) # 6sigma_d


point_group_symmetry_names = ['8C3','8C3','8C3','8C3','8C3','8C3','8C3','8C3',
                              '3C2','3C2','3C2',
                              '6S4','6S4','6S4','6S4','6S4','6S4',
                              '6sigma_d','6sigma_d','6sigma_d','6sigma_d','6sigma_d','6sigma_d',
                              ]


Save_And_Read_SymOperations.SaveSymmetryFile('Td',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_Td')[1]


Inp_Atoms,Structure = load_xyz(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\ModelPolyhedron_8-9-10\CN4_Tetrahedron.xyz',2)
test_pointgroup(Operations,Structure,Inp_Atoms)

#%%
# =============================================================================
# =============================================================================
# # pointgroup_Th
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 8C3
# =============================================================================
axis = [1,1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

# =============================================================================
# 3C2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]])))

# =============================================================================
# 8S6
# =============================================================================
axis = [1,1,1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6
axis = [1,1,1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,1,-1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,1,-1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,-1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

axis = [1,-1,-1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6

# =============================================================================
# 3sigma_h
# =============================================================================

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h

axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h
axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h


point_group_symmetry_names = ['8C3','8C3','8C3','8C3','8C3','8C3','8C3','8C3',
                              '3C2','3C2','3C2',
                              'i',
                              '8S6','8S6','8S6','8S6','8S6','8S6','8S6','8S6',
                              '3sigma_h','3sigma_h','3sigma_h',
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('Th',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_Th')[1]


#%%
# =============================================================================
# =============================================================================
# # pointgroup_S6
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# C3(z)
# =============================================================================
axis = [0,0,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C3(z)
# =============================================================================
# C3(z)^2
# =============================================================================
axis = [0,0,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C3(z)^2

# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]])))
# =============================================================================
# S6
# =============================================================================
axis = [0,0,1]
angle = 60
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6
# =============================================================================
# S6^2
# =============================================================================
axis = [0,0,1]
angle = 60 * 5
Symmetry_Operations.append(Improper_rotation(axis, angle)) # 8S6


point_group_symmetry_names = ['8C3','8C3','8C3','8C3','8C3','8C3','8C3','8C3',
                              '3C2','3C2','3C2',
                              'i',
                              '8S6','8S6','8S6','8S6','8S6','8S6','8S6','8S6',
                              '3sigma_h','3sigma_h','3sigma_h',
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('S6',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_S6')[1]

#%%
# =============================================================================
# =============================================================================
# # pointgroup_D2h
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 3C2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]])))
# =============================================================================
# 3sigma_h
# =============================================================================

axis =[1,0,0]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h
axis =[0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h
axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h


point_group_symmetry_names = ['3C2','3C2','3C2',
                              'i',
                              '3sigma_h','3sigma_h','3sigma_h',
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('D2h',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D2h')[1]
#%%
# =============================================================================
# =============================================================================
# # pointgroup_C2h
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 3C2
# =============================================================================
axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2
# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]])))
# =============================================================================
# 3sigma_h
# =============================================================================

axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h


point_group_symmetry_names = ['C2',
                              'i',
                              'sigma_h',
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('C2h',point_group_symmetry_names,Symmetry_Operations)
symnames, Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C2h')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_S4
# =============================================================================
# =============================================================================
Symmetry_Operations = []

axis = [0,0,1]
angle = 90*2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # C2
angle = 90
Symmetry_Operations.append(Improper_rotation(axis, angle)) # S4
angle = 90 * 3
Symmetry_Operations.append(Improper_rotation(axis, angle)) # S4^3


point_group_symmetry_names = [
    'S4',
    'S4^3',                          
    'C2']
Save_And_Read_SymOperations.SaveSymmetryFile('S4',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_S4')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_T
# =============================================================================
# =============================================================================
Symmetry_Operations = []

# =============================================================================
# 8C3
# =============================================================================
axis = [1,1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

axis = [1,-1,-1]
angle = 120 * 2
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 8C3

# =============================================================================
# 3C2
# =============================================================================
axis = [1,0,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

axis = [0,1,0]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2

axis = [0,0,1]
angle = 180
Symmetry_Operations.append(rotation_matrix(axis, angle)) # 3C2


point_group_symmetry_names = ['8C3','8C3','8C3','8C3','8C3','8C3','8C3','8C3',
                              '3C2','3C2','3C2',

                              ]


Save_And_Read_SymOperations.SaveSymmetryFile('T',point_group_symmetry_names,Symmetry_Operations)
Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_T')[1]


Inp_Atoms,Structure = load_xyz(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\ModelPolyhedron_8-9-10\CN4_Tetrahedron.xyz',2)
test_pointgroup(Operations,Structure,Inp_Atoms)

#%%
# =============================================================================
# =============================================================================
# # pointgroup_D2
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 1C2(z)

axis = [0,1,0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 1C2(y)

axis = [1,0,0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 1C2(x)

point_group_symmetry_names = ['1C2(z)', '1C2(y)', '1C2(x)']
Save_And_Read_SymOperations.SaveSymmetryFile('D2',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D2')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_D5h
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 5 
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
angle = 360 / 5 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5

axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2
angle = 360 / 5 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2


axis = [1,0,0] 
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.cos(2*np.pi / 5),np.sin(2*np.pi / 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.cos(2*np.pi * 2/5),np.sin(2*np.pi * 2/5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.cos(2*np.pi * 3/5),np.sin(2*np.pi * 3/5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.cos(2*np.pi * 4/5),np.sin(2*np.pi * 4/5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2

axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h


axis = [0,0,1]
angle = 360 / 5 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S5
angle = 360 / 5 * 4
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S5

axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S5^3
angle = 360 / 5 *3
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S5^3


axis = [0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis = [np.sin(2*np.pi / 5),np.cos(2*np.pi / 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis = [np.sin(2*np.pi * 2/5),np.cos(2*np.pi * 2/5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis = [np.sin(2*np.pi * 3/5),np.cos(2*np.pi * 3/5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v
axis = [np.sin(2*np.pi * 4/5),np.cos(2*np.pi * 4/5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_v



point_group_symmetry_names = ['2C5', '2C5', '2C5^2','2C5^2',
                              '5C2','5C2','5C2','5C2','5C2',
                              'sigma_h',
                              '2S5','2S5','2S5^3','2S5^3',
                              '5sigma_v','5sigma_v','5sigma_v','5sigma_v','5sigma_v']
Save_And_Read_SymOperations.SaveSymmetryFile('D5h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D5h')

#%%

# =============================================================================
# =============================================================================
# # pointgroup_D5d
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 5 
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
angle = 360 / 5 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2
angle = 360 / 5 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2


axis = [0,1,0] 
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi / 5),np.cos(2*np.pi / 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *2 / 5),np.cos(2*np.pi *2/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *3 / 5),np.cos(2*np.pi *3/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *4 / 5),np.cos(2*np.pi *4/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2

Symmetry_Operations.append(np.array(np.array([[-1, 0, 0], 
                                       [0, -1, 0],
                                       [0, 0, -1]])))   #i

axis = [0,0,1]
angle = 360 * 3 / 10 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10^3
angle = 360 * 7 / 10
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10^3
axis = [0,0,1]
angle = 360 * 5 / 10 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10
angle = 360 * 1 / 10 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10

axis = [0,1,0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_d
axis = [np.sin(2*np.pi / 5),np.cos(2*np.pi / 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_d
axis = [np.sin(2*np.pi *2 / 5),np.cos(2*np.pi *2/ 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_d
axis = [np.sin(2*np.pi *3 / 5),np.cos(2*np.pi *3/ 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_d
axis = [np.sin(2*np.pi *4 / 5),np.cos(2*np.pi *4/ 5),0]
Symmetry_Operations.append(mirror_image(axis)) # 5sigma_d

point_group_symmetry_names = ['2C5', '2C5', '2C5^2','2C5^2',
                              '5C2','5C2','5C2','5C2','5C2',
                              'i',
                              ' 2S10^3',' 2S10^3','2S10','2S10',
                              '5sigma_d','5sigma_d','5sigma_d','5sigma_d','5sigma_d']
Save_And_Read_SymOperations.SaveSymmetryFile('D5d',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D5d')
#%%
# =============================================================================
# =============================================================================
# # pointgroup_D7h
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 7
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C7
angle = 360 / 7 * 6
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C7
angle = 360 / 7 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C7^2
angle = 360 / 7 * 5
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C7^2
angle = 360 / 7 *3 
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C7^3
angle = 360 / 7 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C7^3




axis = [1,0,0] 
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2
axis = [np.cos(2*np.pi / 7),np.sin(2*np.pi / 7),0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2
axis = [np.cos(2*np.pi / 7 * 2),np.sin(2*np.pi / 7* 2),0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2
axis = [np.cos(2*np.pi / 7 * 3),np.sin(2*np.pi / 7* 3),0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2
axis = [np.cos(2*np.pi / 7 * 4),np.sin(2*np.pi / 7* 4),0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2
axis = [np.cos(2*np.pi / 7 * 5),np.sin(2*np.pi / 7* 5),0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2
axis = [np.cos(2*np.pi / 7* 6),np.sin(2*np.pi / 7* 6),0]
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 7C'2



axis =[0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h


axis = [0,0,1]
angle = 360 / 7
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2S7
angle = 360 / 7 * 6
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2S7
angle = 360 / 7 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2S7^3
angle = 360 / 7 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2S7^3
angle = 360 / 7 * 5
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2S7^5
angle = 360 / 7 * 2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2S7^5





axis = [0,1,0] 
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v
axis = [np.sin(2*np.pi / 7),np.cos(2*np.pi / 7),0]
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v
axis = [np.sin(2*np.pi / 7 * 2),np.cos(2*np.pi / 7* 2),0]
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v
axis = [np.sin(2*np.pi / 7 * 3),np.cos(2*np.pi / 7* 3),0]
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v
axis = [np.sin(2*np.pi / 7 * 4),np.cos(2*np.pi / 7* 4),0]
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v
axis = [np.sin(2*np.pi / 7 * 5),np.cos(2*np.pi / 7* 5),0]
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v
axis = [np.sin(2*np.pi / 7* 6),np.cos(2*np.pi / 7* 6),0]
Symmetry_Operations.append( mirror_image(axis)) # 7sigma_v




point_group_symmetry_names = ['2C7', '2C7', 
                              '2C7^2','2C7^2',
                              '2C7^3','2C7^3',
                              '7C_2','7C_2','7C_2','7C_2','7C_2','7C_2','7C_2',
                              'sigma_h',
                              '2S7','2S7',
                              '2S7^3','2S7^3',
                              '2S7^5','2S7^5',
                              '7sigma_v','7sigma_v','7sigma_v','7sigma_v','7sigma_v','7sigma_v','7sigma_v']

Save_And_Read_SymOperations.SaveSymmetryFile('D7h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D7h')
test_pointgroup(Operations,Structure,Inp_Atoms,dev_name)
dev_name, Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D7h')
Inp_Atoms,Structure = load_xyz(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\ModelPolyhedron_8-9-10\CN9_HeptagonalBiPyramid.xyz',2)
test_pointgroup(Operations,Structure,Inp_Atoms,dev_name)
#%%

# =============================================================================
# =============================================================================
# # pointgroup_S10
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 5 
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
angle = 360 / 5 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2
angle = 360 / 5 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2


axis = [0,1,0] 
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi / 5),np.cos(2*np.pi / 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *2 / 5),np.cos(2*np.pi *2/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *3 / 5),np.cos(2*np.pi *3/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *4 / 5),np.cos(2*np.pi *4/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2

Symmetry_Operations.append(np.array(np.array([[-1, 0, 0], 
                                       [0, -1, 0],
                                       [0, 0, -1]])))   #i

axis = [0,0,1]
angle = 360 * 3 / 10 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10^3
angle = 360 * 7 / 10
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10^3
axis = [0,0,1]
angle = 360 * 5 / 10 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10
angle = 360 * 1 / 10 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S10


point_group_symmetry_names = ['2C5', '2C5', '2C5^2','2C5^2',
                              '5C2','5C2','5C2','5C2','5C2',
                              'i',
                              ' 2S10^3',' 2S10^3','2S10','2S10']
Save_And_Read_SymOperations.SaveSymmetryFile('S10',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_S10')
#%%

# =============================================================================
# =============================================================================
# # pointgroup_D5
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 5 
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
angle = 360 / 5 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5
axis = [0,0,1]
angle = 360 / 5 *2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2
angle = 360 / 5 *3
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C5^2


axis = [0,1,0] 
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi / 5),np.cos(2*np.pi / 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *2 / 5),np.cos(2*np.pi *2/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *3 / 5),np.cos(2*np.pi *3/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2
axis = [np.sin(2*np.pi *4 / 5),np.cos(2*np.pi *4/ 5),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 5C'2


point_group_symmetry_names = ['2C5', '2C5', '2C5^2','2C5^2',
                              '5C2','5C2','5C2','5C2','5C2',
                              'i'
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('D5',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D5')


#%%

# =============================================================================
# =============================================================================
# # pointgroup_D8h
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 360 / 8
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C8
angle = 360 / 8 * 7
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C8

axis = [0,0,1]
angle = 360 / 8 * 3 
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C8^3
angle = 360 / 8 * 5
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C8^3

axis = [0,0,1]
angle = 360 / 8 * 2
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C4
angle = 360 / 8 * 6
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 2C4

axis = [0,0,1]
angle = 360 / 8 * 4
Symmetry_Operations.append( rotation_matrix(axis, angle)) # C2

axis = [0,1,0] 
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C_2
axis = [np.sin(2*np.pi / 8),np.cos(2*np.pi / 8),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C_2
axis = [np.sin(2*np.pi *2/ 8),np.cos(2*np.pi *2/ 8),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C_2
axis = [np.sin(2*np.pi *3/ 8),np.cos(2*np.pi *3/ 8),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C_2


axis = [np.sin(2*np.pi / 16),np.cos(2*np.pi / 16),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C__2
axis = [np.sin(2*np.pi *3/ 16),np.cos(2*np.pi *3/ 16),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C__2
axis = [np.sin(2*np.pi *5/ 16),np.cos(2*np.pi *5/ 16),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C__2
axis = [np.sin(2*np.pi *7/ 16),np.cos(2*np.pi *7/ 16),0]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 4C__2


Symmetry_Operations.append(np.array(np.array([[-1, 0, 0], 
                                       [0, -1, 0],
                                       [0, 0, -1]])))   #i

axis = [0,0,1]
angle = 360 / 8
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S8^3
angle = 360 / 8 * 7
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S8^3

axis = [0,0,1]
angle = 360 / 8 * 3 
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S8
angle = 360 / 8 * 5
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S8

axis = [0,0,1]
angle = 360 / 8 * 2
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S4
angle = 360 / 8 * 6
Symmetry_Operations.append( Improper_rotation(axis, angle)) # 2S4


axis = [0,0,1]
Symmetry_Operations.append(mirror_image(axis)) # sigma_h

axis = [0,1,0] 
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_v
axis = [np.sin(2*np.pi / 8),np.cos(2*np.pi / 8),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_v
axis = [np.sin(2*np.pi *2/ 8),np.cos(2*np.pi *2/ 8),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_v
axis = [np.sin(2*np.pi *3/ 8),np.cos(2*np.pi *3/ 8),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_v


axis = [np.sin(2*np.pi / 16),np.cos(2*np.pi / 16),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_d
axis = [np.sin(2*np.pi *3/ 16),np.cos(2*np.pi *3/ 16),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_d
axis = [np.sin(2*np.pi *5/ 16),np.cos(2*np.pi *5/ 16),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_d
axis = [np.sin(2*np.pi *7/ 16),np.cos(2*np.pi *7/ 16),0]
Symmetry_Operations.append( mirror_image(axis)) # 4sigma_d




point_group_symmetry_names = ['2C8', '2C8', '2C8^3','2C8^3',
                              '2C4','2C4',
                              'C2',
                              '4C_2','4C_2','4C_2','4C_2',
                              '4C__2','4C__2','4C__2','4C__2',
                              'i',
                              '2S8^3','2S8^3','2S8','2S8',
                              '2S4','2S4',
                              'sigma_h',
                              'sigma_v','sigma_v','sigma_v','sigma_v',
                              'sigma_d','sigma_d','sigma_d','sigma_d',
                              ]
Save_And_Read_SymOperations.SaveSymmetryFile('D8h',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_D8h')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_C2v
# =============================================================================
# =============================================================================
Symmetry_Operations = []
axis = [0,0,1]
angle = 180
Symmetry_Operations.append( rotation_matrix(axis, angle)) # 1C2(z)

axis = [0,1,0]
Symmetry_Operations.append( mirror_image(axis)) # sigmav(xz)

axis = [1,0,0]
Symmetry_Operations.append( mirror_image(axis)) # sigmav(yz)

point_group_symmetry_names = ['1C2(z)', 'sigmav(xz)', 'sigmav(yz)']
Save_And_Read_SymOperations.SaveSymmetryFile('C2v',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_C2v')


#%%
# =============================================================================
# =============================================================================
# # pointgroup_Ci
# =============================================================================
# =============================================================================


Symmetry_Operations = []
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0], 
                                       [0, -1, 0],
                                       [0, 0, -1]])))   #i


point_group_symmetry_names = ['i']
Save_And_Read_SymOperations.SaveSymmetryFile('Ci',point_group_symmetry_names,Symmetry_Operations)
Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_Ci')

Symmetry_Operations = []


#%%
# =============================================================================
# =============================================================================
# # pointgroup_Ih
# =============================================================================
# =============================================================================
# =============================================================================
# 12C5
# =============================================================================
Symmetry_Operations = []
angle = 360 / 5
X = (1+np.sqrt(5)) / 2
Axis_12_0 = []
Axis_12_0.append([0,1,X])
Axis_12_0.append([0,1,-X])
Axis_12_0.append([0,-1,X])
Axis_12_0.append([0,-1,-X])
Axis_12_0.append([1,X,0])
Axis_12_0.append([1,-X,0])
Axis_12_0.append([-1,X,0])
Axis_12_0.append([-1,-X,0])
Axis_12_0.append([X,0,1])
Axis_12_0.append([X,0,-1])
Axis_12_0.append([-X,0,1])
Axis_12_0.append([-X,0,-1])

Axis_12 = []
R_to_zaxis = rotation_matrix_from_vectors(Axis_12_0[0],[0,0,1])
for axis in Axis_12_0:
    Axis_12.append(np.dot(axis,R_to_zaxis))

for axis in Axis_12:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle)) # 20C5
for axis in Axis_12:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle * 2)) # 20C3^5
# =============================================================================
# 20C3
# =============================================================================
angle = 120
def find_coordinate(point1, point2, point3):
    x_coordinate = (point1[0] + point2[0] + point3[0]) / 3
    y_coordinate = (point1[1] + point2[1] + point3[1]) / 3
    z_coordinate = (point1[2] + point2[2] + point3[2]) / 3
    return [x_coordinate, y_coordinate, z_coordinate]
Axis_20  = []
Axis_20.append(find_coordinate( Axis_12[0] , Axis_12[4], Axis_12[8]  ) )
Axis_20.append(find_coordinate( Axis_12[7] , Axis_12[10], Axis_12[11]  ) )
Axis_20.append(find_coordinate( Axis_12[10] , Axis_12[2], Axis_12[7]  ) )
Axis_20.append(find_coordinate( Axis_12[5] , Axis_12[9], Axis_12[8]  ) )
Axis_20.append(find_coordinate( Axis_12[8] , Axis_12[2], Axis_12[0]  ) )
Axis_20.append(find_coordinate( Axis_12[10] , Axis_12[0], Axis_12[2]  ) )
Axis_20.append(find_coordinate( Axis_12[10] , Axis_12[0], Axis_12[6]  ) )
Axis_20.append(find_coordinate( Axis_12[7] , Axis_12[2], Axis_12[5]  ) )
Axis_20.append(find_coordinate( Axis_12[8] , Axis_12[2], Axis_12[5]  ) )
Axis_20.append(find_coordinate( Axis_12[4] , Axis_12[0], Axis_12[6]  ) )
for axis in Axis_20:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle)) # 20C3
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle*2)) # 20C3
# =============================================================================
# 15C2
# =============================================================================
def find_coordinate(point1, point2):
    x_coordinate = (point1[0] + point2[0]) / 2
    y_coordinate = (point1[1] + point2[1]) / 2
    z_coordinate = (point1[2] + point2[2]) / 2
    return [x_coordinate, y_coordinate, z_coordinate]
Axis_15 = []
Axis_15.append(find_coordinate( Axis_12[6] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[6] , Axis_12[4]) )
Axis_15.append(find_coordinate( Axis_12[4] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[5] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[7] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[8] , Axis_12[9]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[6]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[4]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[2]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[5]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[7]) )
angle = 180
for axis in Axis_15:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle)) # 20C2
# =============================================================================
# i
# =============================================================================
Symmetry_Operations.append(np.array(np.array([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]])))
# =============================================================================
# 12S10
# =============================================================================

angle = 360 / 10
for axis in Axis_12:
    Symmetry_Operations.append(Improper_rotation(np.array(axis), angle)) # 12S10
for axis in Axis_12:
    Symmetry_Operations.append(Improper_rotation(-np.array(axis), angle * 3)) # 12S10^3

# =============================================================================
# 20S6
# =============================================================================
angle = 360 / 6
for axis in Axis_20:
    Symmetry_Operations.append(Improper_rotation(np.array(axis), angle)) # 20S6
    Symmetry_Operations.append(Improper_rotation(np.array(axis), angle*5)) # 20S6
# =============================================================================
# 15sigma
# =============================================================================

for axis in Axis_15:
    Symmetry_Operations.append(mirror_image(np.array(axis))) # 15sigma

point_group_symmetry_names = ['12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5',
                              '12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2',
                              '20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3',
                              '15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2',
                              'i',
                              '12S10','12S10','12S10','12S10','12S10','12S10','12S10','12S10','12S10','12S10','12S10','12S10',
                              '12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3','12S10^3',
                              '20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6','20S6',
                              '15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma','15sigma'
                                                           ]


Save_And_Read_SymOperations.SaveSymmetryFile('Ih',point_group_symmetry_names,Symmetry_Operations)
dev_name,Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_Ih')


#%%

atoms = ['O' for i in range(len(Axis_12))] + ['N' for i in range(len(Axis_20))] +['B' for i in range(len(Axis_15))]

Axis = np.append(Axis_12,Axis_20,axis=0) 
Axis = np.append(Axis,Axis_15,axis=0)  * 2.5

new_coords = Axis.astype(str)
new_coords = np.insert(new_coords, 0, atoms , axis=1)


output_file_path = r'C:\Users\Villads\Documents\Project_SoftwareSymmetryDeviation\Structures_in_different_point_groups\pointgroup_Ih_TRIAL.xyz'

np.savetxt(output_file_path, new_coords, delimiter='\t', fmt='%s', header=str(len(new_coords))+'\n', comments='')

#%%
# =============================================================================
# =============================================================================
# # pointgroup_I
# =============================================================================
# =============================================================================
# =============================================================================
# 12C5
# =============================================================================
Symmetry_Operations = []
angle = 360 / 5
X = (1+np.sqrt(5)) / 2
Axis_12_0 = []
Axis_12_0.append([0,1,X])
Axis_12_0.append([0,1,-X])
Axis_12_0.append([0,-1,X])
Axis_12_0.append([0,-1,-X])
Axis_12_0.append([1,X,0])
Axis_12_0.append([1,-X,0])
Axis_12_0.append([-1,X,0])
Axis_12_0.append([-1,-X,0])
Axis_12_0.append([X,0,1])
Axis_12_0.append([X,0,-1])
Axis_12_0.append([-X,0,1])
Axis_12_0.append([-X,0,-1])

Axis_12 = []
R_to_zaxis = rotation_matrix_from_vectors(Axis_12_0[0],[0,0,1])
for axis in Axis_12_0:
    Axis_12.append(np.dot(axis,R_to_zaxis))

for axis in Axis_12:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle)) # 20C5
for axis in Axis_12:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle * 2)) # 20C3^5
# =============================================================================
# 20C3
# =============================================================================
angle = 120
def find_coordinate(point1, point2, point3):
    x_coordinate = (point1[0] + point2[0] + point3[0]) / 3
    y_coordinate = (point1[1] + point2[1] + point3[1]) / 3
    z_coordinate = (point1[2] + point2[2] + point3[2]) / 3
    return [x_coordinate, y_coordinate, z_coordinate]
Axis_20  = []
Axis_20.append(find_coordinate( Axis_12[0] , Axis_12[4], Axis_12[8]  ) )
Axis_20.append(find_coordinate( Axis_12[7] , Axis_12[10], Axis_12[11]  ) )
Axis_20.append(find_coordinate( Axis_12[10] , Axis_12[2], Axis_12[7]  ) )
Axis_20.append(find_coordinate( Axis_12[5] , Axis_12[9], Axis_12[8]  ) )
Axis_20.append(find_coordinate( Axis_12[8] , Axis_12[2], Axis_12[0]  ) )
Axis_20.append(find_coordinate( Axis_12[10] , Axis_12[0], Axis_12[2]  ) )
Axis_20.append(find_coordinate( Axis_12[10] , Axis_12[0], Axis_12[6]  ) )
Axis_20.append(find_coordinate( Axis_12[7] , Axis_12[2], Axis_12[5]  ) )
Axis_20.append(find_coordinate( Axis_12[8] , Axis_12[2], Axis_12[5]  ) )
Axis_20.append(find_coordinate( Axis_12[4] , Axis_12[0], Axis_12[6]  ) )
for axis in Axis_20:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle)) # 20C3
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle*2)) # 20C3
# =============================================================================
# 15C2
# =============================================================================
def find_coordinate(point1, point2):
    x_coordinate = (point1[0] + point2[0]) / 2
    y_coordinate = (point1[1] + point2[1]) / 2
    z_coordinate = (point1[2] + point2[2]) / 2
    return [x_coordinate, y_coordinate, z_coordinate]
Axis_15 = []
Axis_15.append(find_coordinate( Axis_12[6] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[6] , Axis_12[4]) )
Axis_15.append(find_coordinate( Axis_12[4] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[5] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[7] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[8] , Axis_12[9]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[6]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[4]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[2]) )
Axis_15.append(find_coordinate( Axis_12[0] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[10]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[8]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[5]) )
Axis_15.append(find_coordinate( Axis_12[2] , Axis_12[7]) )
angle = 180
for axis in Axis_15:
    Symmetry_Operations.append(rotation_matrix(np.array(axis), angle)) # 20C2

point_group_symmetry_names = ['12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5','12C5',
                              '12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2','12C5^2',
                              '20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3','20C3',
                              '15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2','15C2',
                                                           ]


Save_And_Read_SymOperations.SaveSymmetryFile('I',point_group_symmetry_names,Symmetry_Operations)
dev_name,Operations = Save_And_Read_SymOperations.ReadSymmetryFile('pointgroup_I')




