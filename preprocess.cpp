#include <fstream>
#include <cstring>
#include <string>
#include <iostream>

#include "common.h"

using namespace std;

int main() {
    cout << "Main start" << endl;
    
    ifstream seq_name_file{"seq_path.txt", ios::in};

    ofstream length_file{"length.data", ios::out};
    FILE* data_file = fopen("seq.data", "w");

    string path_str;
    char seq_genome[MAX_SEQUENCE_LENGTH];

    cout << "Start loading" << endl;
    while (getline(seq_name_file, path_str)) {
        ifstream fasta_file{path_str, ios::in};
        cout << "Start process " << path_str << endl;
        string content;

        unsigned int length = 0;

        // 消耗掉首行
        getline(fasta_file, content);

        memset(seq_genome, 0,  MAX_SEQUENCE_LENGTH * sizeof(char));

        while (getline(fasta_file, content)) {
            if (content.length() <= 1) {
                break;
            }

            memcpy(seq_genome + length, content.c_str(), content.length());
            length += content.length();
        }
        length_file << length << endl;
        fwrite(seq_genome,  MAX_SEQUENCE_LENGTH * sizeof(char), 1, data_file);
        fasta_file.close();
        cout << path_str << "\tDone" << endl;
    }

    length_file.close();
    seq_name_file.close();

    fclose(data_file);
    return 0;
}