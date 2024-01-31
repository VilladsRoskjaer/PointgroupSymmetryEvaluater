import csv
import numpy as np
import os

def SaveSymmetryFile(PointGroup,OperationNames,Operations):
    # Saving the lists to a CSV file
    data = list(zip(Operations,OperationNames))
    filename = 'pointgroup_{}.csv'.format(PointGroup)
    save_folderpath = save_folderpath = os.path.join(r'C:\Users\Villads\Documents\Project_SymmetryDeviation\SymmetryOperationLibrary', 'Operations')
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
    
    Path_to_folder = r'C:\Users\Villads\Documents\Project_SymmetryDeviation\SymmetryOperationLibrary\Operations'
    
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
