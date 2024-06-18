
#include "hdf_post_process.h"

#define FORCE_FLOAT_DATA

void process_hdf_file(std::string input_filename, std::string output_filename)
{
    try
    {
        // Open the input file with read-only access
        H5::H5File infile(input_filename, H5F_ACC_RDONLY);

        // Create the output file with H5F_ACC_TRUNC access (this erases all data in an existing file or creates the file if it doesn't exist)
        H5::H5File outfile(output_filename, H5F_ACC_TRUNC);
        H5::Group out_group = outfile.createGroup("/Table_Data");

        // Find all of the channel names
        // Also get the trigger timestamps from the group for the first enabled channel
        // TODO may need to handle this differently when multiple devices exist
        std::map<std::string, hid_t> instr_grp_id_map;
        std::map<std::string, hid_t> chan_grp_id_map;
        std::map<std::string, hid_t> dio_grp_id_map;
        //std::map<std::string, std::vector<std::string>> instr_chan_map;
        H5::Group root = infile.openGroup(root_group_name);
        hsize_t num_instr_objs = root.getNumObjs();
        for(int i = 0; i < num_instr_objs; i++)
        {
            std::string instr_obj_name = root.getObjnameByIdx(i);
            instr_grp_id_map[instr_obj_name] = root.getObjId(instr_obj_name);
            H5O_type_t instr_obj_type = root.childObjType(instr_obj_name);
            switch(instr_obj_type)
            {
            case H5O_TYPE_GROUP:
            {
                //instr_chan_map[instr_obj_name] = std::vector<std::string>();
                H5::Group instrument = root.openGroup(instr_obj_name);
                hsize_t num_chan_objs = instrument.getNumObjs();
                for(int j = 0; j < num_chan_objs; j++)
                {
                    std::string ch_obj_name = instrument.getObjnameByIdx(j);
                    if(ch_obj_name.find("CH") != ch_obj_name.npos)
                    {
                        //instr_chan_map[instr_obj_name].push_back(ch_obj_name);
                        std::string full_ch_name = get_full_channel_name(instr_obj_name, ch_obj_name);
                        chan_grp_id_map[full_ch_name] = instrument.getObjId(ch_obj_name);
                        if(j == 0)
                        {
                            H5::Group ch_grp = instrument.openGroup(ch_obj_name);
                            H5::DataSet context_dset = ch_grp.openDataSet(ex_context_name);
                            write_trigger_timestamps(context_dset, out_group);
                        }
                    }
                    else if(ch_obj_name.find("DIO") != ch_obj_name.npos)
                    {
                        std::string full_ch_name = get_full_channel_name(instr_obj_name, ch_obj_name);
                        dio_grp_id_map[full_ch_name] = instrument.getObjId(ch_obj_name);
                    }
                }
                break;
            }
            default:
                break;
            }
        }

#ifdef FORCE_FLOAT_DATA
        std::map<std::string, std::vector<RecordDataFloat>> chan_data = read_data_float(chan_grp_id_map);
        write_record_table_float(out_group, instr_grp_id_map, chan_grp_id_map, chan_data);
#else
        //std::map<std::string, std::vector<RecordData>> chan_data = read_data(root, instr_chan_map);
        std::map<std::string, std::vector<RecordData>> chan_data = read_data(chan_grp_id_map);
        write_record_table(out_group, instr_grp_id_map, chan_grp_id_map, chan_data);
#endif

        if(!dio_grp_id_map.empty())
        {
            // If there's DIO data, write that out as well
            std::map<std::string, std::vector<DioData>> data = read_dio_data(dio_grp_id_map);
            write_dio_table(out_group, instr_grp_id_map, dio_grp_id_map, data);
        }

        outfile.close();
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

void write_trigger_timestamps(H5::DataSet dset, H5::Group output_grp)
{
    try
    {
        // Get the size of "ContextFields" in the file's compound type
        // We know that the field is a 1-dimensional array (from previously inspecting the type),
        // but we need to know the length for when we process the data after reading it
        H5::CompType datatype_in = dset.getCompType();
        int context_index = datatype_in.getMemberIndex(context_words_field);
        H5::ArrayType arrtype = datatype_in.getMemberArrayType(context_index);
        int arr_ndims = arrtype.getArrayNDims();
        hsize_t* arr_dims = new hsize_t[arr_ndims];
        arrtype.getArrayDims(arr_dims);
        // If the field were a multi-dimensional array, we would need the total size
        size_t arr_len = 1;
        for(int i = 0; i < arr_ndims; i++)
        {
            arr_len *= (size_t)arr_dims[i];
        }

        // Construct the HDF datatype for reading the dataset into memory
        // The member names need to match the fields in the file's compound type
        // Only those members will be retrieved when we call the dataset's read method
        H5::CompType memtype(sizeof(ex_meas_info_dataset_t));
        memtype.insertMember(time_seconds_field, HOFFSET(ex_meas_info_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_upper_field, HOFFSET(ex_meas_info_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_lower_field, HOFFSET(ex_meas_info_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
        memtype.insertMember(context_indicator_field, HOFFSET(ex_meas_info_dataset_t, context_indicator_field), H5::PredType::NATIVE_INT32);
        memtype.insertMember(context_words_field, HOFFSET(ex_meas_info_dataset_t, context_fields), H5::ArrayType(H5::PredType::NATIVE_INT32, arr_ndims, arr_dims));

        // Determine the dataspace dimensions so that we can allocate enough space for all of the data
        H5::DataSpace dspace_in = dset.getSpace();
        int ndims_in = dspace_in.getSimpleExtentNdims();
        hsize_t* dims_in = new hsize_t[ndims_in];
        dspace_in.getSimpleExtentDims(dims_in);

        size_t total_data_size = 1;
        for(int i = 0; i < ndims_in; i++)
        {
            total_data_size *= (size_t)dims_in[i];
        }

        // Read the data
        std::vector<ex_meas_info_dataset_t> dset_vec;
        dset_vec.resize(total_data_size);
        dset.read(dset_vec.data(), memtype);

        // Process the data a little
        std::vector<ExMeasInfo> context_vec = convert_file_data(dset_vec);

        // Create a buffer for writing the output
        std::vector<double> data_out;
        for(size_t i = 0; i < context_vec.size(); i++)
        {
            // Only record trigger timestamps (ignore packets that are CJC values)
            ExMeasInfo pkt = context_vec[i];
            if(pkt.context_indicators.ext_meas_info.trigger)
            {
                double val = context_vec[i].trigger_seconds + context_vec[i].trigger_fraction;
                data_out.emplace_back(val);
            }
        }

        // Create the dataset for writing the output
        hsize_t dims_out[1] = { data_out.size() };
        H5::DataSpace dspace_out = H5::DataSpace(1, dims_out);
        // TODO set DSetCreatPropList to allow for compression?
        H5::DataSet ds_out = output_grp.createDataSet("Trigger Timestamps", H5::PredType::NATIVE_DOUBLE, dspace_out);
        ds_out.write(data_out.data(), H5::PredType::NATIVE_DOUBLE);
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void write_record_table(H5::Group output_grp, std::map<std::string, hid_t> instr_grp_id_map,
    std::map<std::string, hid_t> chan_grp_id_map, std::map<std::string, std::vector<RecordData>> chan_data_map)
{
    try
    {
        size_t num_columns = chan_data_map.size() + 1;
        size_t num_rows = 0;

        std::vector<std::string> sorted_chans = sorted_channel_names(chan_grp_id_map);

        // Use the first channel to determine the relative timestamps of each sample, the total number of samples, and other acquisition info
        std::vector<RecordData> ch_recs = chan_data_map[sorted_chans.at(0)];
        RecordData first_rec = ch_recs.at(0);
        double initial_time = first_rec.time_seconds + first_rec.time_fraction;

        // Record attributes for the group
        int ch_count = (int)chan_grp_id_map.size();
        H5::Attribute ch_count_attr = output_grp.createAttribute(attr_name_chan_count, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        ch_count_attr.write(H5::PredType::NATIVE_INT, &ch_count);
        ch_count_attr.close();

        int num_recs = (int)ch_recs.size();
        H5::Attribute rec_count_attr = output_grp.createAttribute(attr_name_rec_count, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        rec_count_attr.write(H5::PredType::NATIVE_INT, &num_recs);
        rec_count_attr.close();

        int record_size = (int)first_rec.samples.size();
        H5::Attribute rec_size_attr = output_grp.createAttribute(attr_name_rec_size, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        rec_size_attr.write(H5::PredType::NATIVE_INT, &record_size);
        rec_size_attr.close();

        double sample_interval;
        double sample_rate_val;
        H5::Group first_chan = H5::Group((*chan_grp_id_map.begin()).second);
        if(first_chan.attrExists(attr_name_sample_rate))
        {
            H5::Attribute ch_rate_attr = first_chan.openAttribute(attr_name_sample_rate);
            ch_rate_attr.read(H5::PredType::NATIVE_DOUBLE, &sample_rate_val);
            ch_rate_attr.close();

            sample_interval = 1 / sample_rate_val;
        }
        else
        {
            // If the Sample Rate attribute doesn't exist, we can infer the correct value from the record timestamps
            RecordData second_rec = ch_recs.at(1);
            double record_interval = (second_rec.time_seconds - first_rec.time_seconds) + (second_rec.time_fraction - first_rec.time_fraction);
            sample_interval = record_interval / record_size;
            sample_rate_val = 1.0 / sample_interval;
        }
        H5::Attribute sample_rate_attr = output_grp.createAttribute(attr_name_sample_rate + " (Hz)", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
        sample_rate_attr.write(H5::PredType::NATIVE_DOUBLE, &sample_rate_val);
        sample_rate_attr.close();

        num_rows = num_recs * record_size;
        std::vector<double> times = calculate_sample_timestamps(ch_recs, sample_interval);

        // Create the dataset
        hsize_t dims_out[RECORD_TABLE_NUM_DIMS] = { num_rows, num_columns };
        hsize_t cdims[RECORD_TABLE_NUM_DIMS] = { num_rows, num_columns };
        H5::DataSpace dspace_out = H5::DataSpace(RECORD_TABLE_NUM_DIMS, dims_out);
        // use a dataset creation prop list to enable compression
        H5::DSetCreatPropList ds_creatplist;
        ds_creatplist.setChunk(RECORD_TABLE_NUM_DIMS, cdims);
        ds_creatplist.setDeflate(5);
        H5::DataSet ds_out = output_grp.createDataSet(dataset_name_records, H5::PredType::NATIVE_DOUBLE, dspace_out, ds_creatplist);

        // Record attributes for the dataset
        H5::StrType time_col_attr_type = H5::StrType(0, time_col_attr_val.length());
        H5::Attribute time_col_attr = ds_out.createAttribute("Column 0", time_col_attr_type, H5::DataSpace(H5S_SCALAR));
        time_col_attr.write(time_col_attr_type, time_col_attr_val);
        time_col_attr.close();

        H5::Attribute interval_attr = ds_out.createAttribute(attr_name_interval, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
        interval_attr.write(H5::PredType::NATIVE_DOUBLE, &sample_interval);
        interval_attr.close();

        std::string systime_value = "Unknown";
        if(!instr_grp_id_map.empty())
        {
            hid_t first_instr_id = (*instr_grp_id_map.begin()).second;
            H5::Group first_instr = H5::Group(first_instr_id);
            if(first_instr.attrExists(attr_name_system_time))
            {
                H5::Attribute systime_attr = first_instr.openAttribute(attr_name_system_time);
                H5::DataType systime_attr_type = systime_attr.getDataType();
                systime_attr.read(systime_attr_type, systime_value);
                systime_attr.close();
            }
        }
        H5::StrType systime_source_attr_type = H5::StrType(0, systime_value.length());
        H5::Attribute systime_source_attr = ds_out.createAttribute(attr_name_system_time, systime_source_attr_type, H5::DataSpace(H5S_SCALAR));
        systime_source_attr.write(systime_source_attr_type, systime_value);
        systime_source_attr.close();

        if(systime_value == "PTP")
        {
            // When the system time source isn't based on PTP, we can't convert the initial time to a wall clock time
            std::time_t time_val = (std::time_t)(initial_time - LEAP_SECONDS);
            char time_buf[100];
            if(std::strftime(time_buf, sizeof(time_buf), TIME_FORMAT, std::gmtime(&time_val)))
            {
                std::string time_str(time_buf);
                H5::StrType wall_time_attr_type = H5::StrType(0, time_str.length());
                H5::Attribute wall_time_attr = ds_out.createAttribute(attr_name_initial_time_wall, wall_time_attr_type, H5::DataSpace(H5S_SCALAR));
                wall_time_attr.write(wall_time_attr_type, time_str);
                wall_time_attr.close();
            }
        }
        // We always want to record the initial time in seconds (even when we have the wall clock time as well)
        H5::Attribute time_zero_attr = ds_out.createAttribute(attr_name_initial_time_secs, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
        time_zero_attr.write(H5::PredType::NATIVE_DOUBLE, &initial_time);
        time_zero_attr.close();

        // Define the memory dataspace and hyperslab to write the sample interval times
        hsize_t col_dims[1] = { times.size() };
        H5::DataSpace times_dspace(1, col_dims);
        hsize_t count[RECORD_TABLE_NUM_DIMS] = { times.size(), 1 };
        hsize_t offset[RECORD_TABLE_NUM_DIMS] = { 0, 0 };
        dspace_out.selectHyperslab(H5S_SELECT_SET, count, offset);
        ds_out.write(times.data(), H5::PredType::NATIVE_DOUBLE, times_dspace, dspace_out);

        // Write the channel data
        // Each channel's data will be the same length as the sample times column
        hsize_t ch_col = 1;
        for(auto it = sorted_chans.begin(); it != sorted_chans.end(); it++)
        {
            // Record channel attributes
            std::stringstream ss;
            ss << "Column " << ch_col;
            std::string ch_col_attr_name = ss.str();

            std::string ch_name = *it;
            H5::Group ch_group = H5::Group(chan_grp_id_map[ch_name]);
            std::string ch_desc = get_channel_description(ch_group);

            std::string ch_col_description = ch_name + " (" + ch_desc + ")";
            H5::StrType ch_col_attr_type = H5::StrType(0, ch_col_description.length());
            H5::Attribute ch_col_attr = ds_out.createAttribute(ch_col_attr_name, ch_col_attr_type, H5::DataSpace(H5S_SCALAR));
            ch_col_attr.write(ch_col_attr_type, ch_col_description);
            ch_col_attr.close();

            // Write out the channel data
#ifdef WRITE_SINGLE_RECORD
            // TODO FIXME this only works for the first record (not sure why)
            H5::DataSpace ch_dspace(1, col_dims);
            write_channel_data(ds_out, dspace_out, ch_dspace, ch_col, chan_data_map[ch_name];
#else
            write_channel_data(ds_out, dspace_out, ch_col, chan_data_map[ch_name]);
#endif
            ch_col++;
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

void write_channel_data(H5::DataSet dset_out, H5::DataSpace dspace_out, hsize_t ch_col_index, std::vector<RecordData> ch_data)
{
    try
    {
        // Flatten the channel data
        std::vector<double> samples;
        for(auto i = ch_data.begin(); i != ch_data.end(); i++)
        {
            for(auto j = (*i).samples.begin(); j != (*i).samples.end(); j++)
            {
                samples.emplace_back(*j);
            }
        }

        // Write the channel data
        hsize_t col_dims[1] = { samples.size() };
        H5::DataSpace dspace_ch(1, col_dims);

        hsize_t count[RECORD_TABLE_NUM_DIMS] = { samples.size(), 1 };
        hsize_t offset[RECORD_TABLE_NUM_DIMS] = { 0, ch_col_index };
        dspace_out.selectHyperslab(H5S_SELECT_SET, count, offset);
        dset_out.write(samples.data(), H5::PredType::NATIVE_DOUBLE, dspace_ch, dspace_out);
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

void write_channel_data(H5::DataSet dset_out, H5::DataSpace dspace_out, H5::DataSpace dspace_ch, hsize_t ch_col_index, std::vector<RecordData> ch_data)
{
    try
    {
        // Write the channel data one record at a time
        // TODO FIXME this only works for the first record (not sure why)
        hsize_t row_offset = 0;
        for(size_t i = 0; i < ch_data.size(); i++)
        {
            RecordData rec = ch_data[i];

            hsize_t count[RECORD_TABLE_NUM_DIMS] = { rec.samples.size(), 1 };
            hsize_t offset[RECORD_TABLE_NUM_DIMS] = { row_offset, ch_col_index };
            dspace_out.selectHyperslab(H5S_SELECT_SET, count, offset);

            hsize_t ch_count[1] = { rec.samples.size() };
            hsize_t ch_offset[1] = { row_offset };
            dspace_ch.selectHyperslab(H5S_SELECT_SET, ch_count, ch_offset);

            dset_out.write(rec.samples.data(), H5::PredType::NATIVE_DOUBLE, dspace_ch, dspace_out);
            row_offset += rec.samples.size();
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

std::vector<RecordData> read_channel_data(H5::DataSet dset, std::string dset_name)
{
    std::vector<RecordData> retval;
    try
    {
        // Get the size of the "Samples" field in the file's compound type
        H5::CompType datatype = dset.getCompType();
        int samp_index = datatype.getMemberIndex(samples_field);
        H5::ArrayType arrtype = datatype.getMemberArrayType(samp_index);
        int arr_ndims = arrtype.getArrayNDims();
        hsize_t* arr_dims = new hsize_t[arr_ndims];
        arrtype.getArrayDims(arr_dims);
        // If the field were a multi-dimensional array, we would need the total size
        size_t arr_len = 1;
        for(int i = 0; i < arr_ndims; i++)
        {
            arr_len *= (size_t)arr_dims[i];
        }

        // Determine the dataspace dimensions so that we can allocate enough space for all of the data
        H5::DataSpace dspace = dset.getSpace();
        int ndims = dspace.getSimpleExtentNdims();
        hsize_t* dims = new hsize_t[ndims];
        dspace.getSimpleExtentDims(dims);

        size_t total_data_size = 1;
        for(int i = 0; i < ndims; i++)
        {
            total_data_size *= (size_t)dims[i];
        }

        // Construct the HDF datatype for reading the dataset into memory
        if(dset_name == meas_float_name)
        {
            H5::CompType memtype(sizeof(meas_float_dataset_t));
            memtype.insertMember(time_seconds_field, HOFFSET(meas_float_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_upper_field, HOFFSET(meas_float_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_lower_field, HOFFSET(meas_float_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
            memtype.insertMember(samples_field, HOFFSET(meas_float_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_FLOAT, arr_ndims, arr_dims));

            std::vector<meas_float_dataset_t> data_vec;
            data_vec.resize(total_data_size);
            dset.read(data_vec.data(), memtype);
            retval = convert_file_data(data_vec, arr_len);
        }
        else if(dset_name == meas_float_name)
        {
            H5::CompType memtype(sizeof(meas_fixed_point_dataset_t));
            memtype.insertMember(time_seconds_field, HOFFSET(meas_fixed_point_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_upper_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_lower_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
            memtype.insertMember(samples_field, HOFFSET(meas_fixed_point_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_INT32, arr_ndims, arr_dims));

            std::vector<meas_fixed_point_dataset_t> data_vec;
            data_vec.resize(total_data_size);
            dset.read(data_vec.data(), memtype);
            retval = convert_file_data(data_vec, arr_len);
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

std::map<std::string, std::vector<RecordData>> read_data(H5::Group root, std::map<std::string, std::vector<std::string>> instr_chan_map)
{
    std::map<std::string, std::vector<RecordData>> retval;
    try
    {
        for(auto it = instr_chan_map.begin(); it != instr_chan_map.end(); it++)
        {
            std::string instr_name = (*it).first;
            std::vector<std::string> instr_chans = (*it).second;
            H5::Group instrument = root.openGroup(instr_name);
            for(size_t i = 0; i < instr_chans.size(); i++)
            {
                std::string ch = instr_chans[i];
                std::string full_ch_name = get_full_channel_name(instr_name, ch);
                H5::Group ch_grp = instrument.openGroup(ch);
                hsize_t num_ch_objs = ch_grp.getNumObjs();
                for(int j = 0; j < num_ch_objs; j++)
                {
                    std::string ds_name = ch_grp.getObjnameByIdx(j);
                    // skip the data sets containing context packets
                    if(ds_name.find("INFO") == ds_name.npos)
                    {
                        H5::DataSet ds = ch_grp.openDataSet(ds_name);
                        retval[full_ch_name] = read_channel_data(ds, ds_name);
                    }
                }
            }
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

std::map<std::string, std::vector<RecordData>> read_data(std::map<std::string, hid_t> chan_grp_id_map)
{
    std::map<std::string, std::vector<RecordData>> retval;
    try
    {
        for(auto it = chan_grp_id_map.begin(); it != chan_grp_id_map.end(); it++)
        {
            std::string full_ch_name = (*it).first;
            H5::Group ch_grp = H5::Group((*it).second);
            hsize_t num_ch_objs = ch_grp.getNumObjs();
            for(int j = 0; j < num_ch_objs; j++)
            {
                std::string ds_name = ch_grp.getObjnameByIdx(j);
                // skip the data sets containing context packets
                if(ds_name.find("INFO") == ds_name.npos)
                {
                    H5::DataSet ds = ch_grp.openDataSet(ds_name);
                    retval[full_ch_name] = read_channel_data(ds, ds_name);
                }
            }
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

void write_record_table_float(H5::Group output_grp, std::map<std::string, hid_t> instr_grp_id_map,
    std::map<std::string, hid_t> chan_grp_id_map, std::map<std::string, std::vector<RecordDataFloat>> chan_data_map)
{
    try
    {
        size_t num_columns = chan_data_map.size() + 1;
        size_t num_rows = 0;

        std::vector<std::string> sorted_chans = sorted_channel_names(chan_grp_id_map);

        // Use the first channel to determine the relative timestamps of each sample, the total number of samples, and other acquisition info
        std::vector<RecordDataFloat> ch_recs = chan_data_map[sorted_chans.at(0)];
        RecordDataFloat first_rec = ch_recs.at(0);
        double initial_time = first_rec.time_seconds + first_rec.time_fraction;

        // Record attributes for the group
        int ch_count = (int)chan_grp_id_map.size();
        H5::Attribute ch_count_attr = output_grp.createAttribute(attr_name_chan_count, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        ch_count_attr.write(H5::PredType::NATIVE_INT, &ch_count);
        ch_count_attr.close();

        int num_recs = (int)ch_recs.size();
        H5::Attribute rec_count_attr = output_grp.createAttribute(attr_name_rec_count, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        rec_count_attr.write(H5::PredType::NATIVE_INT, &num_recs);
        rec_count_attr.close();

        int record_size = (int)first_rec.samples.size();
        H5::Attribute rec_size_attr = output_grp.createAttribute(attr_name_rec_size, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        rec_size_attr.write(H5::PredType::NATIVE_INT, &record_size);
        rec_size_attr.close();

        double sample_interval;
        double sample_rate_val;
        H5::Group first_chan = H5::Group((*chan_grp_id_map.begin()).second);
        if(first_chan.attrExists(attr_name_sample_rate))
        {
            H5::Attribute ch_rate_attr = first_chan.openAttribute(attr_name_sample_rate);
            ch_rate_attr.read(H5::PredType::NATIVE_DOUBLE, &sample_rate_val);
            ch_rate_attr.close();

            sample_interval = 1 / sample_rate_val;
        }
        else
        {
            // If the Sample Rate attribute doesn't exist, we can infer the correct value from the record timestamps
            RecordDataFloat second_rec = ch_recs.at(1);
            double record_interval = (second_rec.time_seconds - first_rec.time_seconds) + (second_rec.time_fraction - first_rec.time_fraction);
            sample_interval = record_interval / record_size;
            sample_rate_val = 1.0 / sample_interval;
        }
        H5::Attribute sample_rate_attr = output_grp.createAttribute(attr_name_sample_rate + " (Hz)", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
        sample_rate_attr.write(H5::PredType::NATIVE_DOUBLE, &sample_rate_val);
        sample_rate_attr.close();

        num_rows = num_recs * record_size;
        std::vector<float> times = calculate_sample_timestamps_float(ch_recs, sample_interval);

        // Create the dataset
        hsize_t dims_out[RECORD_TABLE_NUM_DIMS] = { num_rows, num_columns };
        hsize_t cdims[RECORD_TABLE_NUM_DIMS] = { num_rows, num_columns };
        H5::DataSpace dspace_out = H5::DataSpace(RECORD_TABLE_NUM_DIMS, dims_out);
        // use a dataset creation prop list to enable compression
        H5::DSetCreatPropList ds_creatplist;
        ds_creatplist.setChunk(RECORD_TABLE_NUM_DIMS, cdims);
        ds_creatplist.setDeflate(5);
        H5::DataSet ds_out = output_grp.createDataSet(dataset_name_records, H5::PredType::NATIVE_FLOAT, dspace_out, ds_creatplist);

        // Record attributes for the dataset
        H5::StrType time_col_attr_type = H5::StrType(0, time_col_attr_val.length());
        H5::Attribute time_col_attr = ds_out.createAttribute("Column 0", time_col_attr_type, H5::DataSpace(H5S_SCALAR));
        time_col_attr.write(time_col_attr_type, time_col_attr_val);
        time_col_attr.close();

        H5::Attribute interval_attr = ds_out.createAttribute(attr_name_interval, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
        interval_attr.write(H5::PredType::NATIVE_DOUBLE, &sample_interval);
        interval_attr.close();

        std::string systime_value = "Unknown";
        if(!instr_grp_id_map.empty())
        {
            hid_t first_instr_id = (*instr_grp_id_map.begin()).second;
            H5::Group first_instr = H5::Group(first_instr_id);
            if(first_instr.attrExists(attr_name_system_time))
            {
                H5::Attribute systime_attr = first_instr.openAttribute(attr_name_system_time);
                H5::DataType systime_attr_type = systime_attr.getDataType();
                systime_attr.read(systime_attr_type, systime_value);
                systime_attr.close();
            }
        }
        H5::StrType systime_source_attr_type = H5::StrType(0, systime_value.length());
        H5::Attribute systime_source_attr = ds_out.createAttribute(attr_name_system_time, systime_source_attr_type, H5::DataSpace(H5S_SCALAR));
        systime_source_attr.write(systime_source_attr_type, systime_value);
        systime_source_attr.close();

        if(systime_value == "PTP")
        {
            // When the system time source isn't based on PTP, we can't convert the initial time to a wall clock time
            std::time_t time_val = (std::time_t)(initial_time - LEAP_SECONDS);
            char time_buf[100];
            if(std::strftime(time_buf, sizeof(time_buf), TIME_FORMAT, std::gmtime(&time_val)))
            {
                std::string time_str(time_buf);
                H5::StrType wall_time_attr_type = H5::StrType(0, time_str.length());
                H5::Attribute wall_time_attr = ds_out.createAttribute(attr_name_initial_time_wall, wall_time_attr_type, H5::DataSpace(H5S_SCALAR));
                wall_time_attr.write(wall_time_attr_type, time_str);
                wall_time_attr.close();
            }
        }
        // We always want to record the initial time in seconds (even when we have the wall clock time as well)
        H5::Attribute time_zero_attr = ds_out.createAttribute(attr_name_initial_time_secs, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
        time_zero_attr.write(H5::PredType::NATIVE_DOUBLE, &initial_time);
        time_zero_attr.close();

        // Define the memory dataspace and hyperslab to write the sample interval times
        hsize_t col_dims[1] = { times.size() };
        H5::DataSpace times_dspace(1, col_dims);
        hsize_t count[RECORD_TABLE_NUM_DIMS] = { times.size(), 1 };
        hsize_t offset[RECORD_TABLE_NUM_DIMS] = { 0, 0 };
        dspace_out.selectHyperslab(H5S_SELECT_SET, count, offset);
        ds_out.write(times.data(), H5::PredType::NATIVE_FLOAT, times_dspace, dspace_out);

        // Write the channel data
        // Each channel's data will be the same length as the sample times column
        hsize_t ch_col = 1;
        for(auto it = sorted_chans.begin(); it != sorted_chans.end(); it++)
        {
            // Record channel attributes
            std::stringstream ss;
            ss << "Column " << ch_col;
            std::string ch_col_attr_name = ss.str();

            std::string ch_name = *it;
            H5::Group ch_group = H5::Group(chan_grp_id_map[ch_name]);
            std::string ch_desc = get_channel_description(ch_group);

            std::string ch_col_description = ch_name + " (" + ch_desc + ")";
            H5::StrType ch_col_attr_type = H5::StrType(0, ch_col_description.length());
            H5::Attribute ch_col_attr = ds_out.createAttribute(ch_col_attr_name, ch_col_attr_type, H5::DataSpace(H5S_SCALAR));
            ch_col_attr.write(ch_col_attr_type, ch_col_description);
            ch_col_attr.close();

            // Write out the channel data
            write_channel_data_float(ds_out, dspace_out, ch_col, chan_data_map[ch_name]);
            ch_col++;
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

void write_channel_data_float(H5::DataSet dset_out, H5::DataSpace dspace_out, hsize_t ch_col_index, std::vector<RecordDataFloat> ch_data)
{
    try
    {
        // Flatten the channel data
        std::vector<float> samples;
        for(auto i = ch_data.begin(); i != ch_data.end(); i++)
        {
            for(auto j = (*i).samples.begin(); j != (*i).samples.end(); j++)
            {
                samples.emplace_back(*j);
            }
        }

        // Write the channel data
        hsize_t col_dims[1] = { samples.size() };
        H5::DataSpace dspace_ch(1, col_dims);

        hsize_t count[RECORD_TABLE_NUM_DIMS] = { samples.size(), 1 };
        hsize_t offset[RECORD_TABLE_NUM_DIMS] = { 0, ch_col_index };
        dspace_out.selectHyperslab(H5S_SELECT_SET, count, offset);
        dset_out.write(samples.data(), H5::PredType::NATIVE_FLOAT, dspace_ch, dspace_out);
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

std::vector<RecordDataFloat> read_channel_data_float(H5::DataSet dset, std::string dset_name)
{
    std::vector<RecordDataFloat> retval;
    try
    {
        // Get the size of the "Samples" field in the file's compound type
        H5::CompType datatype = dset.getCompType();
        int samp_index = datatype.getMemberIndex(samples_field);
        H5::ArrayType arrtype = datatype.getMemberArrayType(samp_index);
        int arr_ndims = arrtype.getArrayNDims();
        hsize_t* arr_dims = new hsize_t[arr_ndims];
        arrtype.getArrayDims(arr_dims);
        // If the field were a multi-dimensional array, we would need the total size
        size_t arr_len = 1;
        for(int i = 0; i < arr_ndims; i++)
        {
            arr_len *= (size_t)arr_dims[i];
        }

        // Determine the dataspace dimensions so that we can allocate enough space for all of the data
        H5::DataSpace dspace = dset.getSpace();
        int ndims = dspace.getSimpleExtentNdims();
        hsize_t* dims = new hsize_t[ndims];
        dspace.getSimpleExtentDims(dims);

        size_t total_data_size = 1;
        for(int i = 0; i < ndims; i++)
        {
            total_data_size *= (size_t)dims[i];
        }

        // Construct the HDF datatype for reading the dataset into memory
        H5::CompType memtype(sizeof(meas_float_dataset_t));
        memtype.insertMember(time_seconds_field, HOFFSET(meas_float_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_upper_field, HOFFSET(meas_float_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_lower_field, HOFFSET(meas_float_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
        memtype.insertMember(samples_field, HOFFSET(meas_float_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_FLOAT, arr_ndims, arr_dims));

        std::vector<meas_float_dataset_t> data_vec;
        data_vec.resize(total_data_size);
        dset.read(data_vec.data(), memtype);
        if(dset_name == meas_float_name)
        {
            H5::CompType memtype(sizeof(meas_float_dataset_t));
            memtype.insertMember(time_seconds_field, HOFFSET(meas_float_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_upper_field, HOFFSET(meas_float_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_lower_field, HOFFSET(meas_float_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
            memtype.insertMember(samples_field, HOFFSET(meas_float_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_FLOAT, arr_ndims, arr_dims));

            std::vector<meas_float_dataset_t> data_vec;
            data_vec.resize(total_data_size);
            dset.read(data_vec.data(), memtype);
            retval = convert_file_data_float(data_vec, arr_len);
        }
        else if(dset_name == meas_float_name)
        {
            H5::CompType memtype(sizeof(meas_fixed_point_dataset_t));
            memtype.insertMember(time_seconds_field, HOFFSET(meas_fixed_point_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_upper_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
            memtype.insertMember(time_pico_lower_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
            memtype.insertMember(samples_field, HOFFSET(meas_fixed_point_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_INT32, arr_ndims, arr_dims));

            std::vector<meas_fixed_point_dataset_t> data_vec;
            data_vec.resize(total_data_size);
            dset.read(data_vec.data(), memtype);
            retval = convert_file_data_float(data_vec, arr_len); // TODO check whether we lose too much precision here
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

std::map<std::string, std::vector<RecordDataFloat>> read_data_float(H5::Group root, std::map<std::string, std::vector<std::string>> instr_chan_map)
{
    std::map<std::string, std::vector<RecordDataFloat>> retval;
    try
    {
        for(auto it = instr_chan_map.begin(); it != instr_chan_map.end(); it++)
        {
            std::string instr_name = (*it).first;
            std::vector<std::string> instr_chans = (*it).second;
            H5::Group instrument = root.openGroup(instr_name);
            for(size_t i = 0; i < instr_chans.size(); i++)
            {
                std::string ch = instr_chans[i];
                std::string full_ch_name = get_full_channel_name(instr_name, ch);
                H5::Group ch_grp = instrument.openGroup(ch);
                hsize_t num_ch_objs = ch_grp.getNumObjs();
                for(int j = 0; j < num_ch_objs; j++)
                {
                    std::string ds_name = ch_grp.getObjnameByIdx(j);
                    // skip the data sets containing context packets
                    if(ds_name.find("INFO") == ds_name.npos)
                    {
                        H5::DataSet ds = ch_grp.openDataSet(ds_name);
                        retval[full_ch_name] = read_channel_data_float(ds, ds_name);
                    }
                }
            }
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

std::map<std::string, std::vector<RecordDataFloat>> read_data_float(std::map<std::string, hid_t> chan_grp_id_map)
{
    std::map<std::string, std::vector<RecordDataFloat>> retval;
    try
    {
        for(auto it = chan_grp_id_map.begin(); it != chan_grp_id_map.end(); it++)
        {
            std::string full_ch_name = (*it).first;
            H5::Group ch_grp = H5::Group((*it).second);
            hsize_t num_ch_objs = ch_grp.getNumObjs();
            for(int j = 0; j < num_ch_objs; j++)
            {
                std::string ds_name = ch_grp.getObjnameByIdx(j);
                // skip the data sets containing context packets
                if(ds_name.find("INFO") == ds_name.npos)
                {
                    H5::DataSet ds = ch_grp.openDataSet(ds_name);
                    retval[full_ch_name] = read_channel_data_float(ds, ds_name);
                }
            }
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

std::vector<DioData> read_dio_data(H5::DataSet dset)
{
    std::vector<DioData> retval;
    try
    {
        // Get the size of the "Samples" field in the file's compound type
        H5::CompType datatype = dset.getCompType();
        int samp_index = datatype.getMemberIndex(samples_field);
        H5::ArrayType arrtype = datatype.getMemberArrayType(samp_index);
        int arr_ndims = arrtype.getArrayNDims();
        hsize_t* arr_dims = new hsize_t[arr_ndims];
        arrtype.getArrayDims(arr_dims);
        // If the field were a multi-dimensional array, we would need the total size
        size_t arr_len = 1;
        for(int i = 0; i < arr_ndims; i++)
        {
            arr_len *= (size_t)arr_dims[i];
        }

        // Determine the dataspace dimensions so that we can allocate enough space for all of the data
        H5::DataSpace dspace = dset.getSpace();
        int ndims = dspace.getSimpleExtentNdims();
        hsize_t* dims = new hsize_t[ndims];
        dspace.getSimpleExtentDims(dims);

        size_t total_data_size = 1;
        for(int i = 0; i < ndims; i++)
        {
            total_data_size *= (size_t)dims[i];
        }

        // Construct the HDF datatype for reading the dataset into memory
        H5::CompType memtype(sizeof(meas_fixed_point_dataset_t));
        memtype.insertMember(time_seconds_field, HOFFSET(meas_fixed_point_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_upper_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_lower_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
        memtype.insertMember(samples_field, HOFFSET(meas_fixed_point_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_INT32, arr_ndims, arr_dims));

        std::vector<meas_fixed_point_dataset_t> data_vec;
        data_vec.resize(total_data_size);
        dset.read(data_vec.data(), memtype);
        retval = convert_dio_data(data_vec, arr_len);
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

std::map<std::string, std::vector<DioData>> read_dio_data(std::map<std::string, hid_t> dio_grp_id_map)
{
    std::map<std::string, std::vector<DioData>> retval;
    try
    {
        for(auto it = dio_grp_id_map.begin(); it != dio_grp_id_map.end(); it++)
        {
            std::string full_ch_name = (*it).first;
            H5::Group ch_grp = H5::Group((*it).second);
            hsize_t num_ch_objs = ch_grp.getNumObjs();
            for(int j = 0; j < num_ch_objs; j++)
            {
                std::string ds_name = ch_grp.getObjnameByIdx(j);
                // skip the data sets containing context packets
                if(ds_name.find("INFO") == ds_name.npos)
                {
                    H5::DataSet ds = ch_grp.openDataSet(ds_name);
                    retval[full_ch_name] = read_dio_data(ds);
                }
            }
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
    return retval;
}

void write_dio_table(H5::Group output_grp, std::map<std::string, hid_t> instr_grp_id_map,
    std::map<std::string, hid_t> dio_grp_id_map, std::map<std::string, std::vector<DioData>> dio_data_map)
{
    try
    {
        size_t num_columns = dio_data_map.size();
        size_t num_rows = 0;

        // Use the first channel to determine the relative timestamps of each sample, the total number of samples, and other acquisition info
        std::vector<DioData> dio_recs = (*dio_data_map.begin()).second;
        DioData first_rec = dio_recs.at(0);

        // Record attributes for the group
        int dio_count = (int)dio_data_map.size();
        H5::Attribute dio_count_attr = output_grp.createAttribute(attr_name_dio_count, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        dio_count_attr.write(H5::PredType::NATIVE_INT, &dio_count);
        dio_count_attr.close();

        int num_recs = (int)dio_recs.size();
        int record_size = (int)first_rec.samples.size();
        num_rows = num_recs * record_size;

        // Create the dataset
        hsize_t dims_out[RECORD_TABLE_NUM_DIMS] = { num_rows, num_columns };
        hsize_t cdims[RECORD_TABLE_NUM_DIMS] = { num_rows, num_columns };
        H5::DataSpace dspace_out = H5::DataSpace(RECORD_TABLE_NUM_DIMS, dims_out);
        // use a dataset creation prop list to enable compression
        H5::DSetCreatPropList ds_creatplist;
        ds_creatplist.setChunk(RECORD_TABLE_NUM_DIMS, cdims);
        ds_creatplist.setDeflate(5);
        H5::DataSet ds_out = output_grp.createDataSet(dataset_name_dio, H5::PredType::NATIVE_INT, dspace_out, ds_creatplist);

        // Define the memory dataspace and hyperslab to write the sample interval times
        hsize_t count[RECORD_TABLE_NUM_DIMS] = { num_rows, 1 };
        hsize_t offset[RECORD_TABLE_NUM_DIMS] = { 0, 0 };

        // Write the DIO data
        hsize_t ch_col = 0;
        for(auto it = dio_data_map.begin(); it != dio_data_map.end(); it++)
        {
            // Record channel attributes
            std::stringstream ss;
            ss << "Column " << ch_col;
            std::string ch_col_attr_name = ss.str();

            std::string ch_name = (*it).first;
            H5::Group ch_group = H5::Group(dio_grp_id_map[ch_name]);

            std::string ch_col_description = ch_name;
            H5::StrType ch_col_attr_type = H5::StrType(0, ch_col_description.length());
            H5::Attribute ch_col_attr = ds_out.createAttribute(ch_col_attr_name, ch_col_attr_type, H5::DataSpace(H5S_SCALAR));
            ch_col_attr.write(ch_col_attr_type, ch_col_description);
            ch_col_attr.close();

            // Write out the DIO data
            write_dio_data(ds_out, dspace_out, ch_col, (*it).second);
            ch_col++;
        }
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

void write_dio_data(H5::DataSet dset_out, H5::DataSpace dspace_out, hsize_t dio_col_index, std::vector<DioData> dio_data)
{
    try
    {
        // Flatten the data
        std::vector<int> samples;
        for(auto i = dio_data.begin(); i != dio_data.end(); i++)
        {
            for(auto j = (*i).samples.begin(); j != (*i).samples.end(); j++)
            {
                samples.emplace_back(*j);
            }
        }

        // Write the channel data
        hsize_t col_dims[1] = { samples.size() };
        H5::DataSpace dspace_ch(1, col_dims);

        hsize_t count[RECORD_TABLE_NUM_DIMS] = { samples.size(), 1 };
        hsize_t offset[RECORD_TABLE_NUM_DIMS] = { 0, dio_col_index };
        dspace_out.selectHyperslab(H5S_SELECT_SET, count, offset);
        dset_out.write(samples.data(), H5::PredType::NATIVE_INT, dspace_ch, dspace_out);
    }
    catch(H5::Exception)
    {
        H5::Exception::printErrorStack();
    }
}

std::vector<double> calculate_sample_timestamps(std::vector<RecordData> recs, double sample_interval)
{
    int initial_secs = recs.at(0).time_seconds;
    double initial_frac = recs.at(0).time_fraction;
    std::vector<double> ts;
    for(auto it = recs.begin(); it != recs.end(); it++)
    {
        int relative_secs = (*it).time_seconds - initial_secs;
        double relative_frac = (*it).time_fraction - initial_frac;
        for(size_t j = 0; j < (*it).samples.size(); j++)
        {
            //double ts_val = (relative_secs + relative_frac) + j * sample_interval;
            //ts.emplace_back(ts_val);
            ts.push_back((relative_secs + relative_frac) + j * sample_interval);
        }
    }
    return ts;
}

std::vector<float> calculate_sample_timestamps_float(std::vector<RecordDataFloat> recs, double sample_interval)
{
    int initial_secs = recs.at(0).time_seconds;
    double initial_frac = recs.at(0).time_fraction;
    std::vector<float> ts;
    for(auto it = recs.begin(); it != recs.end(); it++)
    {
        int relative_secs = (*it).time_seconds - initial_secs;
        double relative_frac = (*it).time_fraction - initial_frac;
        for(size_t j = 0; j < (*it).samples.size(); j++)
        {
            double ts_val = (relative_secs + relative_frac) + j * sample_interval;
            ts.push_back((float)ts_val);
            //ts.push_back((relative_secs + relative_frac) + j * sample_interval);
        }
    }
    return ts;
}

int channel_index(const std::string &ch)
{
  size_t pos = ch.find("CH");
  std::string num_str = ch.substr(pos + 2);
  return std::stoi(num_str) - 1;
}

bool CompareChannelNames(const std::string &left, const std::string &right)
{
  return channel_index(left) < channel_index(right);
}

std::vector<std::string> sorted_channel_names(std::map<std::string, hid_t> chan_grp_id_map)
{
    std::vector<std::string> names;
    for(auto it = chan_grp_id_map.begin(); it != chan_grp_id_map.end(); it++)
    {
        names.insert(std::lower_bound(names.begin(), names.end(), it->first, &CompareChannelNames), it->first);
    }
    return names;
}

std::string get_full_channel_name(std::string instr, std::string channel)
{
    // Each instrument name is "Slot X"; we only want the number portion
    // This lets us name the channels as 2!CH1, etc, so they match the repeated capability names used in the VTEXDigitizer/VTEXDsa drivers
    std::string slot_prefix = "";
    size_t pos = instr.find(" ");
    if(pos != instr.npos)
    {
        slot_prefix = instr.substr(pos + 1) + "!";
    }
    return slot_prefix + channel;
}

std::string get_channel_description(H5::Group ch_grp)
{
    std::string desc = "";

    H5::Attribute function_attr = ch_grp.openAttribute("Function");
    H5::DataType function_attr_type = function_attr.getDataType();
    std::string function_value;
    function_attr.read(function_attr_type, function_value);
    function_attr.close();
    desc.append(function_value);

    if(function_value == "Temperature")
    {
        if(ch_grp.attrExists("Temperature Units"))
        {
            H5::Attribute units_attr = ch_grp.openAttribute("Temperature Units");
            H5::DataType units_attr_type = units_attr.getDataType();
            std::string units_value;
            units_attr.read(units_attr_type, units_value);
            units_attr.close();
            desc.append(", ");
            desc.append(units_value);
        }

        if(ch_grp.attrExists("Temperature Transducer Type"))
        {
            H5::Attribute transducer_attr = ch_grp.openAttribute("Temperature Transducer Type");
            H5::DataType transducer_attr_type = transducer_attr.getDataType();
            std::string transducer_value;
            transducer_attr.read(transducer_attr_type, transducer_value);
            transducer_attr.close();

            if(transducer_value == "Thermocouple" && ch_grp.attrExists("Thermocouple Type"))
            {
                H5::Attribute tc_attr = ch_grp.openAttribute("Thermocouple Type");
                H5::DataType tc_attr_type = tc_attr.getDataType();
                std::string tc_value;
                tc_attr.read(tc_attr_type, tc_value);
                tc_attr.close();
                desc.append(", Thermocouple Type ");
                desc.append(tc_value);
            }
        }
    }

    return desc;
}
