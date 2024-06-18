#pragma once

#include "hdf_utils.h"

// methods for inspecting the contents of an HDF5 file
void inspect_hdf_file(std::string filename);
void inspect_attributes(H5::Group grp);
void inspect_dataset(H5::DataSet dset);
void inspect_compound_type(H5::CompType ctype);
void inspect_group(H5::Group grp, bool read_attrs = true, bool read_children = true);

// methods for reading data in an HDF5 file
void read_instrument_data(std::string filename);
void read_channel_data(H5::Group grp);

void read_float_dataset(H5::DataSet dset);
std::vector<RecordData> convert_file_data(std::vector<meas_float_dataset_t> data_vec, size_t samples_array_len);
std::vector<RecordDataFloat> convert_file_data_float(std::vector<meas_float_dataset_t> data_vec, size_t samples_array_len);

void read_fixed_point_dataset(H5::DataSet dset);
std::vector<RecordData> convert_file_data(std::vector<meas_fixed_point_dataset_t> data_vec, size_t samples_array_len);
std::vector<RecordDataFloat> convert_file_data_float(std::vector<meas_fixed_point_dataset_t> data_vec, size_t samples_array_len);

void read_ex_meas_info(H5::DataSet dset);
std::vector<ExMeasInfo> convert_file_data(std::vector<ex_meas_info_dataset_t> context_vec);

std::vector<DioData> convert_dio_data(std::vector<meas_fixed_point_dataset_t> data_vec, size_t samples_array_len);
