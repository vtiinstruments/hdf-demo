#pragma once

#include "hdf_utils.h"
#include "hdf_inspection.h"

const std::string time_col_attr_val = "Relative Time (secs)";
const std::string attr_name_interval = "Sample Interval (secs)";
const std::string attr_name_initial_time_secs = "Table Time Zero (secs)";
const std::string attr_name_initial_time_wall = "Table Time Zero";
const std::string attr_name_chan_count = "Enabled Channel Count";
const std::string attr_name_dio_count = "Enabled DIO Count";
const std::string attr_name_rec_count = "Record Count";
const std::string attr_name_rec_size = "Record Size";
const std::string attr_name_sample_rate = "Sample Rate";
const std::string attr_name_system_time = "System Time Source";
const std::string dataset_name_records = "Record Data";
const std::string dataset_name_dio = "DIO Data";

// PTP-based timestamps are in TAI, so the leap seconds must be subtracted before converting to UTC
// This is the current value as 31 December 2016, when the most recent leap second was added
const int LEAP_SECONDS = 37;
#define TIME_FORMAT "%F %T UTC"

#define RECORD_TABLE_NUM_DIMS 2

void process_hdf_file(std::string input_filename, std::string output_filename);

void write_trigger_timestamps(H5::DataSet dset, H5::Group output_grp);

void write_record_table(H5::Group output_grp, std::map<std::string, hid_t> instr_grp_id_map,
    std::map<std::string, hid_t> chan_grp_id_map, std::map<std::string, std::vector<RecordData>> chan_data_map);

void write_channel_data(H5::DataSet dset_out, H5::DataSpace dspace_out, hsize_t ch_col_index, std::vector<RecordData> ch_data);
void write_channel_data(H5::DataSet dset_out, H5::DataSpace dspace_out, H5::DataSpace dspace_ch, hsize_t ch_col_index, std::vector<RecordData> ch_data);

std::vector<RecordData> read_channel_data(H5::DataSet dset, std::string dset_name);
std::map<std::string, std::vector<RecordData>> read_data(H5::Group root, std::map<std::string, std::vector<std::string>> instr_chan_map);
std::map<std::string, std::vector<RecordData>> read_data(std::map<std::string, hid_t> chan_grp_id_map);

void write_record_table_float(H5::Group output_grp, std::map<std::string, hid_t> instr_grp_id_map,
    std::map<std::string, hid_t> chan_grp_id_map, std::map<std::string, std::vector<RecordDataFloat>> chan_data_map);

void write_channel_data_float(H5::DataSet dset_out, H5::DataSpace dspace_out, hsize_t ch_col_index, std::vector<RecordDataFloat> ch_data);

std::vector<RecordDataFloat> read_channel_data_float(H5::DataSet dset, std::string dset_name);
std::map<std::string, std::vector<RecordDataFloat>> read_data_float(H5::Group root, std::map<std::string, std::vector<std::string>> instr_chan_map);
std::map<std::string, std::vector<RecordDataFloat>> read_data_float(std::map<std::string, hid_t> chan_grp_id_map);

std::vector<DioData> read_dio_data(H5::DataSet dset);
std::map<std::string, std::vector<DioData>> read_dio_data(std::map<std::string, hid_t> dio_grp_id_map);

void write_dio_table(H5::Group output_grp, std::map<std::string, hid_t> instr_grp_id_map,
    std::map<std::string, hid_t> dio_grp_id_map, std::map<std::string, std::vector<DioData>> dio_data_map);
void write_dio_data(H5::DataSet dset_out, H5::DataSpace dspace_out, hsize_t dio_col_index, std::vector<DioData> dio_data);

std::vector<double> calculate_sample_timestamps(std::vector<RecordData> recs, double sample_interval);
std::vector<float> calculate_sample_timestamps_float(std::vector<RecordDataFloat> recs, double sample_interval);

std::vector<std::string> sorted_channel_names(std::map<std::string, hid_t> chan_grp_id_map);
std::string get_full_channel_name(std::string instr, std::string channel);
std::string get_channel_description(H5::Group ch_grp);
