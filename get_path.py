import os 
from os import path

dir_path = "G:\\FTP\\TransTemp\\Seme6\\genome_data"

if __name__ == '__main__':
    fasta_file = os.listdir(dir_path)

    with open("seq_path.txt", 'w') as path_file:
        with open("seq_name.txt", 'w') as name_file:
            for file_name in fasta_file:
                if len(file_name) > 6 and file_name[-6:] == '.fasta':
                    path_file.write(path.join(dir_path, file_name) + '\n')
                    name_file.write(file_name[:-6] + '\n')


