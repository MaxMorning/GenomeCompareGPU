#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include "common.h"

#define TILE_SIZE 32

using namespace std;

__constant__ unsigned int compare_seq_length;
__constant__ char compare_seq_genome[MAX_SEQUENCE_LENGTH];

// 计算两条DNA的 Levenshtein 距离
// 通常 str_1 为待比对的DNA，str_2为reference sequence
__device__ int LevenshteinDistance(unsigned int str_2_length, char* str_2)
{
    int buffer_line_0[MAX_SEQUENCE_LENGTH];
    int buffer_line_1[MAX_SEQUENCE_LENGTH];

    int* buffer_line_first = buffer_line_0;
    int* buffer_line_second = buffer_line_1;

    // init first line
    for (int i = 0; i <= compare_seq_length; ++i) {
        buffer_line_first[i] = i;
    }

    for (int line_idx = 1; line_idx <= str_2_length; ++line_idx) {
        buffer_line_second[0] = line_idx;

        for (int i = 1; i <= compare_seq_length; ++i) {
            int replace_delta = 1;
            if (compare_seq_genome[i - 1] == str_2[line_idx - 1]) {
                replace_delta = 0;
            }

            int temp_min = min(min(buffer_line_first[i - 1] + replace_delta, buffer_line_first[i] + 1), buffer_line_second[i - 1] + 1);
            buffer_line_second[i] = temp_min;
        }

        // swap two buffer
        int* temp_buffer_ptr = buffer_line_first;
        buffer_line_first = buffer_line_second;
        buffer_line_second = temp_buffer_ptr;
    }

    int ret_value = buffer_line_first[compare_seq_length];

    return ret_value;
}

// 核函数
// 目前有 1367 条DNA
// 每个block有32 * 32 线程束
__global__ void kernel_func(unsigned int* ref_seq_length_array, char* ref_seq_genome_array, int* distant_result, float* simularity_result) {
    int my_work_seq = threadIdx.x + threadIdx.y * blockDim.x + blockDim.x * blockDim.y * blockIdx.x;
    
    if (my_work_seq >= SEQUENCE_CNT) {
        return;
    }

    char* ref_seq = ref_seq_genome_array + my_work_seq * MAX_SEQUENCE_LENGTH;
    int ref_seq_length = ref_seq_length_array[my_work_seq];
    int distance = LevenshteinDistance(ref_seq_length, ref_seq);
    
    distant_result[my_work_seq] = distance;
    simularity_result[my_work_seq] = 100 - distance * 100.0f / max(ref_seq_length, compare_seq_length);
}

__host__ void load_data_from_disk(unsigned int* h_ref_seq_length_array, char* h_ref_seq_genome_array) {
    ifstream length_file{"length.data", ios::in};
    for (int i = 0; i < SEQUENCE_CNT; ++i) {
        length_file >> h_ref_seq_length_array[i];
    }
    length_file.close();


    FILE* seq_file = fopen("seq.data", "r");
    fread(h_ref_seq_genome_array, MAX_SEQUENCE_LENGTH, SEQUENCE_CNT, seq_file);
    fclose(seq_file);
}

__host__ void load_compare_seq(const string& path) {
    ifstream compare_file{path, ios::in};

    char buffer[MAX_SEQUENCE_LENGTH];
    memset(buffer, 0, MAX_SEQUENCE_LENGTH);

    string buffer_str;
    // 去除首行
    getline(compare_file, buffer_str);

    int char_cnt = 0;
    while (getline(compare_file, buffer_str)) {
        memcpy(buffer + char_cnt, buffer_str.c_str(), buffer_str.length());
        char_cnt += buffer_str.length();
    }

    unsigned int length = strlen(buffer);

    cudaMemcpyToSymbol(compare_seq_length, &length, sizeof(unsigned int));
    cudaMemcpyToSymbol(compare_seq_genome, buffer, length);

    compare_file.close();
}

struct SortableItem {
    int distant;
    float simularity;
    string name;

    bool operator<(const SortableItem& item) const {
        return this->distant < item.distant;
    }
};

__host__ void store_result_to_file(const string& path, int* distant_result, float* simularity_result) {
    ifstream seq_name_file{"seq_name.txt", ios::in};
    ofstream result_file{path, ios::out};
    
    string name;

    vector<SortableItem> result_vector;

    for (int i = 0; i < SEQUENCE_CNT; ++i) {
        getline(seq_name_file, name);
        result_vector.push_back({distant_result[i], simularity_result[i], name});
    }

    sort(result_vector.begin(), result_vector.end());

    for (int i = 0; i < SEQUENCE_CNT; ++i) {
        result_file << setw(16) << result_vector[i].name << setw(8) << result_vector[i].distant << setw(16) << result_vector[i].simularity << endl;
    }

    seq_name_file.close();
    result_file.close();
}

int main() {
    dim3 grid{(SEQUENCE_CNT - 1) / (TILE_SIZE * TILE_SIZE) + 1, 1, 1};
    dim3 block{TILE_SIZE, TILE_SIZE, 1};

    unsigned int* h_ref_seq_length_array = new unsigned int[SEQUENCE_CNT];
    char* h_ref_seq_genome_array = new char[SEQUENCE_CNT * MAX_SEQUENCE_LENGTH];

    load_data_from_disk(h_ref_seq_length_array, h_ref_seq_genome_array);


    unsigned int* d_ref_seq_length_array;
    cudaMalloc(&d_ref_seq_length_array, SEQUENCE_CNT * sizeof(unsigned int));

    char* d_ref_seq_genome_array;
    cudaMalloc(&d_ref_seq_genome_array, SEQUENCE_CNT * MAX_SEQUENCE_LENGTH * sizeof(char));

    cudaMemcpy(d_ref_seq_length_array, h_ref_seq_length_array, SEQUENCE_CNT * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ref_seq_genome_array, h_ref_seq_genome_array, SEQUENCE_CNT * MAX_SEQUENCE_LENGTH * sizeof(char), cudaMemcpyHostToDevice);

    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) 
    {
        fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(cudaStatus));
    }
    else {
        cout << "Copy ref data to device." << endl;
    }

    // 读取待比较序列
    load_compare_seq("sequences.fasta");
    
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) 
    {
        fprintf(stderr, "Load compare sequence failed: %s\n", cudaGetErrorString(cudaStatus));
    }
    else {
        cout << "Load compare data to device." << endl;
    }

    // 为返回值分配显存
    int* d_distant_result;
    cudaMalloc(&d_distant_result, SEQUENCE_CNT * sizeof(int));

    float* d_simularity_result;
    cudaMalloc(&d_simularity_result, SEQUENCE_CNT * sizeof(float));

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) 
    {
        fprintf(stderr, "Alloc result array failed: %s\n", cudaGetErrorString(cudaStatus));
    }
    else {
        cout << "Alloc result array in device." << endl;
    }

    cout << "Start calculating..." << endl;
    kernel_func <<<grid, block>>> (d_ref_seq_length_array, d_ref_seq_genome_array, d_distant_result, d_simularity_result);

    cudaDeviceSynchronize();

    cout << "Calc done." << endl;

    int h_distant_result[SEQUENCE_CNT];
    float h_simularity_result[SEQUENCE_CNT];

    cudaMemcpy(h_distant_result, d_distant_result, SEQUENCE_CNT * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_simularity_result, d_simularity_result, SEQUENCE_CNT * sizeof(float), cudaMemcpyDeviceToHost);
    
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) 
    {
        fprintf(stderr, "copy result to host failed: %s\n", cudaGetErrorString(cudaStatus));
    }
    else {
        cout << "Copy result to host." << endl;
    }
    // 存回磁盘
    store_result_to_file("result.txt", h_distant_result, h_simularity_result);

    cudaFree(d_simularity_result);
    cudaFree(d_distant_result);

    cudaFree(d_ref_seq_genome_array);
    cudaFree(d_ref_seq_length_array);

    delete[] h_ref_seq_genome_array;
    delete[] h_ref_seq_length_array;
    return 0;
}