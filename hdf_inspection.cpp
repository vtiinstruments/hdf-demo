
#include "hdf_inspection.h"

void inspect_hdf_file(std::string filename)
{
    try
    {
        // Open the file with read-only access
        H5::H5File h5file(filename, H5F_ACC_RDONLY);
        std::cout << "HDF5 Filename: " << filename << std::endl;

        // Get the top-level group
        H5::Group root = h5file.openGroup(root_group_name);
        std::cout << "Root Group Name: " << root_group_name << std::endl;

        // Read the top-level group's attributes
        inspect_attributes(root);

        // Read the top-level group's children
        inspect_group(root, false);
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void inspect_attributes(H5::Group grp)
{
    try
    {
        std::cout << "Group Attributes: " << std::endl;
        int group_attrs = grp.getNumAttrs();
        for(int i = 0; i < group_attrs; i++)
        {
            H5::Attribute attr = grp.openAttribute(i);

            // Get the attribute name
            std::string attr_name = attr.getName();

            // Get the generic datatype (in memory)
            H5::DataType attr_type = attr.getDataType();

            // Get the class of the datatype, so we can use an appropriate buffer when calling the read method
            H5T_class_t attr_type_cls = attr.getTypeClass();
            switch(attr_type_cls)
            {
            case H5T_INTEGER:
            {
                int value;
                attr.read(attr_type, &value);
                std::cout << attr_name << " = " << value << std::endl;
                break;
            }
            case H5T_STRING:
            {
                std::string value;
                attr.read(attr_type, value);
                std::cout << attr_name << " = " << value << std::endl;
                break;
            }
            case H5T_FLOAT:
            {
                double value;
                attr.read(attr_type, &value);
                std::cout << attr_name << " = " << value << std::endl;
                break;
            }
            default:
            {
                // Just print the name & type class for other types that require more specific handling
                std::cout << attr_name << ", " << attr_type_cls << std::endl;
            }
            }
        }
        std::cout << std::endl;
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void inspect_dataset(H5::DataSet dset)
{
    try
    {
        // Get the class of the dataset
        H5T_class_t type_class = dset.getTypeClass();

        // Get the dataspace
        H5::DataSpace dspace = dset.getSpace();

        // Get the number of dimensions in the dataspace
        int num_dims = dspace.getSimpleExtentNdims();
        std::cout << "num_dims = " << num_dims << ", ";

        // Get the size of each dimension
        hsize_t* dims = new hsize_t[num_dims];
        int ndims = dspace.getSimpleExtentDims(dims);
        std::cout << "dimensions: ";
        for(int i = 0; i < num_dims; i++)
        {
            if(i > 0) {
                std::cout << " x ";
            }
            std::cout << dims[i];
        }
        std::cout << std::endl;

        // Inspect the contents of the dataspace
        switch(type_class)
        {
        case H5T_COMPOUND:
        {
            H5::CompType mtype = dset.getCompType();
            std::cout << "Class: H5T_COMPOUND" << std::endl;
            inspect_compound_type(mtype);
            break;
        }
        default:
        {
            //TODO
        }
        }

        std::cout << std::endl;
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void inspect_compound_type(H5::CompType ctype)
{
    try
    {
        int nmembs = ctype.getNmembers();
        for(int i = 0; i < nmembs; i++)
        {
            std::string mem_name = ctype.getMemberName(i);
            H5T_class_t mem_class = ctype.getMemberClass(i);
            size_t offset = ctype.getMemberOffset(i);
            std::cout << mem_name << ", type class: " << mem_class << std::endl;
            switch(mem_class)
            {
            case H5T_COMPOUND:
            {
                // Recursively inspect the nested type
                H5::CompType t = ctype.getMemberCompType(i);
                std::cout << std::endl;
                inspect_compound_type(t);
                break;
            }
            case H5T_ARRAY:
            {
                H5::ArrayType t = ctype.getMemberArrayType(i);
                size_t sz = t.getSize();
                H5::DataType base = t.getSuper();
                H5T_class_t base_class = base.getClass();

                std::string base_name = "", ord_str = "";
                hid_t tid = base.getId();
                size_t sz2 = H5Tget_size(tid);
                if(base_class == H5T_INTEGER || base_class == H5T_FLOAT)
                {
                    base_name = base_class == H5T_INTEGER ? "H5T_INTEGER" : "H5T_FLOAT";
                    H5T_order_t ord = H5Tget_order(tid);
                    if(ord == H5T_ORDER_LE)
                        ord_str = "Little endian byte ordering (0)";
                    else if(ord == H5T_ORDER_BE)
                        ord_str = "Big endian byte ordering (1)";
                    else if(ord == H5T_ORDER_VAX)
                        ord_str = "VAX mixed byte ordering (2)";
                }

                printf("  offset: %zd, size: %zd", offset, sz);
                std::cout << ", base type class: " << base_class;
                if(!base_name.empty())
                {
                    std::cout << " (" << base_name << ")";
                }
                if(!ord_str.empty())
                {
                    printf(", base size: %zd, %s", sz2, ord_str.c_str());
                }
                printf("\n");
                break;
            }
            case H5T_INTEGER:
            {
                H5::IntType t = ctype.getMemberIntType(i);
                size_t sz = t.getSize();
                std::string ord_str;
                H5T_order_t ord = t.getOrder(ord_str);
                printf("  offset: %zd, size: %zd, %s\n", offset, sz, ord_str.c_str());
                break;
            }
            case H5T_FLOAT:
            {
                H5::FloatType t = ctype.getMemberFloatType(i);
                size_t sz = t.getSize();
                std::string ord_str;
                H5T_order_t ord = t.getOrder(ord_str);
                printf("  offset: %zd, size: %zd, %s\n", offset, sz, ord_str.c_str());
                break;
            }
            case H5T_STRING:
            {
                // TODO
                break;
            }
            default:
            {
                // TODO
                break;
            }
            }
        }
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void inspect_group(H5::Group grp, bool read_attrs, bool read_children)
{
    try
    {
        if(read_attrs)
        {
            // Read the group's attributes
            inspect_attributes(grp);
        }

        // Read the group's children - these will be either groups or datasets, depending on the file structure
        hsize_t num_objs = grp.getNumObjs();
        std::cout << "Child Objects: " << num_objs << std::endl;
        for(int i = 0; i < num_objs; i++)
        {
            std::string obj_name = grp.getObjnameByIdx(i);
            std::cout << "Object Name: " << obj_name << std::endl;

            if(read_children || obj_name == "CH1")
            {
                H5O_type_t obj_type = grp.childObjType(obj_name);
                switch(obj_type)
                {
                case H5O_TYPE_GROUP:
                {
                    H5::Group child_group = grp.openGroup(obj_name);
                    if(obj_name == "CH1")
                    {
                        inspect_group(child_group, false, true);
                    }
                    else
                    {
                        inspect_group(child_group, true, false);
                    }
                    std::cout << std::endl;
                    break;
                }
                case H5O_TYPE_DATASET:
                {
                    H5::DataSet child_dataset = grp.openDataSet(obj_name);
                    inspect_dataset(child_dataset);
                    break;
                }
                default:
                    // nothing to do here
                    break;
                }
            }
        }
        std::cout << std::endl;
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void read_instrument_data(std::string filename)
{
    try
    {
        // Open the file with read-only access and open the top-level group
        std::cout << "Reading data from: " << filename << std::endl;
        H5::H5File h5file(filename, H5F_ACC_RDONLY);
        H5::Group root = h5file.openGroup(root_group_name);

        // Open the group for each instrument and find the group for CH1
        hsize_t num_objs = root.getNumObjs();
        for(int i = 0; i < num_objs; i++)
        {
            std::string obj_name = root.getObjnameByIdx(i);
            H5O_type_t obj_type = root.childObjType(obj_name);
            switch(obj_type)
            {
            case H5O_TYPE_GROUP:
            {
                H5::Group instrument = root.openGroup(obj_name);
                H5::Group ch1 = instrument.openGroup("CH1");
                read_channel_data(ch1);
                break;
            }
            default:
                // nothing to do here
                break;
            }
        }
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void read_channel_data(H5::Group grp)
{
    try
    {
        // All of the items in this group should be datasets - one for each type of VRT packet received
        hsize_t num_objs = grp.getNumObjs();
        for(int i = 0; i < num_objs; i++)
        {
            std::string obj_name = grp.getObjnameByIdx(i);
            H5O_type_t obj_type = grp.childObjType(obj_name);
            switch(obj_type)
            {
            case H5O_TYPE_DATASET:
            {
                // Datasets with "info" in the name are VRT context packets.
                // Context packets contain auxilary information which may be useful for interpreting the acquisition data.
                if(obj_name == ex_context_name)
                {
                    H5::DataSet channel_dataset = grp.openDataSet(obj_name);
                    read_ex_meas_info(channel_dataset);
                }
                else if(obj_name == meas_float_name)
                {
                    H5::DataSet channel_dataset = grp.openDataSet(obj_name);
                    read_float_dataset(channel_dataset);
                }
                else if(obj_name == meas_fixed_point_name)
                {
                    H5::DataSet channel_dataset = grp.openDataSet(obj_name);
                    read_fixed_point_dataset(channel_dataset);
                }
                break;
            }
            default:
                // nothing to do here
                break;
            }
        }
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

void read_float_dataset(H5::DataSet dset)
{
    try
    {
        // Get the size of the "Samples" field in the file's compound type
        // We know that the field is a 1-dimensional array (from previously inspecting the type),
        // but we need to know the length for when we process the data after reading it
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

        // Construct the HDF datatype for reading the dataset into memory
        // The member names need to match the fields in the file's compound type
        // Only those members will be retrieved when we call the dataset's read method
        H5::CompType memtype(sizeof(meas_float_dataset_t));
        memtype.insertMember(time_seconds_field, HOFFSET(meas_float_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_upper_field, HOFFSET(meas_float_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_lower_field, HOFFSET(meas_float_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
        memtype.insertMember(samples_field, HOFFSET(meas_float_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_FLOAT, arr_ndims, arr_dims));

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

        // Read the data
        std::vector<meas_float_dataset_t> data_vec;
        data_vec.resize(total_data_size);
        dset.read(data_vec.data(), memtype);

        // Process the data a little
        std::vector<RecordData> record_vec = convert_file_data(data_vec, arr_len);

        if(!record_vec.empty())
        {
            printf("Seconds: %d, Fraction: %.9f, Number of Samples: %zd\n", record_vec[0].time_seconds, record_vec[0].time_fraction, record_vec[0].samples.size());
        }
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

inline double convert_meas_int32(int val, const double scale_factor)
{
    return val * scale_factor;
}

inline float convert_meas_int32_float(int val, const double scale_factor)
{
    // TODO FIXME pretty sure that this is a really good way to lose precision
    return (float)(val * scale_factor);
}

std::vector<RecordData> convert_file_data(std::vector<meas_float_dataset_t> data_vec, size_t samples_array_len)
{
    // Convert the data from the format used to store in HDF5 files to something that's a little easier to work with
    // This format is similar to the way data is returned from the Measurement.Read methods in the VTEXDigitizer/VTEXDsa drivers
    std::vector<RecordData> retval;
    for(auto it = data_vec.begin(); it != data_vec.end(); it++)
    {
        RecordData rec;

        // No conversion required for the seconds portion of the timestamp
        rec.time_seconds = (*it).time_sec;

        // Convert the picoseconds fields to fractional seconds
        PICO pico = (PICO)(*it).time_pico_upper; // Don't shift until after we assign to avoid data loss
        pico <<= PICO_SHIFT;
        pico |= (*it).time_pico_lower;
        rec.time_fraction = pico * PICOSECOND_EXPONENT;

        // Add the samples to a vector
        std::vector<double> sample_data;
        sample_data.assign(static_cast<float*>((*it).samples), static_cast<float*>((*it).samples) + samples_array_len);
        rec.samples = sample_data;

        retval.push_back(rec);
    }
    return retval;
}

std::vector<RecordDataFloat> convert_file_data_float(std::vector<meas_float_dataset_t> data_vec, size_t samples_array_len)
{
    // Convert the data from the format used to store in HDF5 files to something that's a little easier to work with
    // This format is similar to the way data is returned from the Measurement.Read methods in the VTEXDigitizer/VTEXDsa drivers
    std::vector<RecordDataFloat> retval;
    for(auto it = data_vec.begin(); it != data_vec.end(); it++)
    {
        RecordDataFloat rec;

        // No conversion required for the seconds portion of the timestamp
        rec.time_seconds = (*it).time_sec;

        // Convert the picoseconds fields to fractional seconds
        PICO pico = (PICO)(*it).time_pico_upper; // Don't shift until after we assign to avoid data loss
        pico <<= PICO_SHIFT;
        pico |= (*it).time_pico_lower;
        rec.time_fraction = pico * PICOSECOND_EXPONENT;

        // Add the samples to a vector
        std::vector<float> sample_data;
        sample_data.assign((*it).samples, (*it).samples + samples_array_len);
        rec.samples = sample_data;

        retval.push_back(rec);
    }
    return retval;
}

void read_fixed_point_dataset(H5::DataSet dset)
{
    try
    {
        // Get the size of the "Samples" field in the file's compound type
        // We know that the field is a 1-dimensional array (from previously inspecting the type),
        // but we need to know the length for when we process the data after reading it
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

        // Construct the HDF datatype for reading the dataset into memory
        // The member names need to match the fields in the file's compound type
        // Only those members will be retrieved when we call the dataset's read method
        H5::CompType memtype(sizeof(meas_fixed_point_dataset_t));
        memtype.insertMember(time_seconds_field, HOFFSET(meas_fixed_point_dataset_t, time_sec), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_upper_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_upper), H5::PredType::NATIVE_INT32);
        memtype.insertMember(time_pico_lower_field, HOFFSET(meas_fixed_point_dataset_t, time_pico_lower), H5::PredType::NATIVE_INT32);
        memtype.insertMember(samples_field, HOFFSET(meas_fixed_point_dataset_t, samples), H5::ArrayType(H5::PredType::NATIVE_INT32, arr_ndims, arr_dims));

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

        // Read the data
        std::vector<meas_fixed_point_dataset_t> data_vec;
        data_vec.resize(total_data_size);
        dset.read(data_vec.data(), memtype);

        // Process the data a little
        std::vector<RecordData> record_vec = convert_file_data(data_vec, arr_len);

        if(!record_vec.empty())
        {
            printf("Seconds: %d, Fraction: %.9f, Number of Samples: %zd\n", record_vec[0].time_seconds, record_vec[0].time_fraction, record_vec[0].samples.size());
        }
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

std::vector<RecordData> convert_file_data(std::vector<meas_fixed_point_dataset_t> data_vec, size_t samples_array_len)
{
    // Convert the data from the format used to store in HDF5 files to something that's a little easier to work with
    // This format is similar to the way data is returned from the Measurement.Read methods in the VTEXDigitizer/VTEXDsa drivers
    std::vector<RecordData> retval;
    double scale_factor = (1.0) / ((double)(1 << MEAS_INT32_SHIFT));
    for(auto it = data_vec.begin(); it != data_vec.end(); it++)
    {
        RecordData rec;

        // No conversion required for the seconds portion of the timestamp
        rec.time_seconds = (*it).time_sec;

        // Convert the picoseconds fields to fractional seconds
        PICO pico = (PICO)(*it).time_pico_upper; // Don't shift until after we assign to avoid data loss
        pico <<= PICO_SHIFT;
        pico |= (*it).time_pico_lower;
        rec.time_fraction = pico * PICOSECOND_EXPONENT;

        // Add the samples to a vector
        std::vector<uint32_t> sample_data;
        sample_data.assign(static_cast<uint32_t*>((*it).samples), static_cast<uint32_t*>((*it).samples) + samples_array_len);

        // if the data format is FixedPoint, it needs to be scaled appropriately
        rec.samples.reserve(sample_data.size());
        for(auto s = sample_data.begin(); s != sample_data.end(); s++)
        {
            rec.samples.push_back(convert_meas_int32((int)*s, scale_factor));
        }

        retval.push_back(rec);
    }
    return retval;
}

std::vector<RecordDataFloat> convert_file_data_float(std::vector<meas_fixed_point_dataset_t> data_vec, size_t samples_array_len)
{
    // Convert the data from the format used to store in HDF5 files to something that's a little easier to work with
    // This format is similar to the way data is returned from the Measurement.Read methods in the VTEXDigitizer/VTEXDsa drivers
    std::vector<RecordDataFloat> retval;
    double scale_factor = (1.0) / ((double)(1 << MEAS_INT32_SHIFT));
    for(auto it = data_vec.begin(); it != data_vec.end(); it++)
    {
        RecordDataFloat rec;

        // No conversion required for the seconds portion of the timestamp
        rec.time_seconds = (*it).time_sec;

        // Convert the picoseconds fields to fractional seconds
        PICO pico = (PICO)(*it).time_pico_upper; // Don't shift until after we assign to avoid data loss
        pico <<= PICO_SHIFT;
        pico |= (*it).time_pico_lower;
        rec.time_fraction = pico * PICOSECOND_EXPONENT;

        // Add the samples to a vector
        std::vector<uint32_t> sample_data;
        sample_data.assign(static_cast<uint32_t*>((*it).samples), static_cast<uint32_t*>((*it).samples) + samples_array_len);

        // if the data format is FixedPoint, it needs to be scaled appropriately
        rec.samples.reserve(sample_data.size());
        for(auto s = sample_data.begin(); s != sample_data.end(); s++)
        {
            // TODO FIXME pretty sure that this is a really good way to lose precision
            rec.samples.push_back(convert_meas_int32_float((int)*s, scale_factor));
        }

        retval.push_back(rec);
    }
    return retval;
}

void read_ex_meas_info(H5::DataSet dset)
{
    try
    {
        // Get the size of "ContextFields" in the file's compound type
        // We know that the field is a 1-dimensional array (from previously inspecting the type),
        // but we need to know the length for when we process the data after reading it
        H5::CompType datatype = dset.getCompType();
        int context_index = datatype.getMemberIndex(context_words_field);
        H5::ArrayType arrtype = datatype.getMemberArrayType(context_index);
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
        H5::DataSpace dspace = dset.getSpace();
        int ndims = dspace.getSimpleExtentNdims();
        hsize_t* dims = new hsize_t[ndims];
        dspace.getSimpleExtentDims(dims);

        size_t total_data_size = 1;
        for(int i = 0; i < ndims; i++)
        {
            total_data_size *= (size_t)dims[i];
        }

        // Read the data
        std::vector<ex_meas_info_dataset_t> dset_vec;
        dset_vec.resize(total_data_size);
        dset.read(dset_vec.data(), memtype);
        
        // Process the data a little
        std::vector<ExMeasInfo> context_vec = convert_file_data(dset_vec);

        if(!context_vec.empty())
        {
            printf("ExMeasInfo, Seconds: %d, Fraction: %.9f\n", context_vec[0].time_seconds, context_vec[0].time_fraction);
        }
    }
    catch(H5::Exception)
    {
        // Print the HDF stack to stderr
        H5::Exception::printErrorStack();
    }
}

std::vector<ExMeasInfo> convert_file_data(std::vector<ex_meas_info_dataset_t> context_vec)
{
    std::vector<ExMeasInfo> retval;
    for(auto it = context_vec.begin(); it != context_vec.end(); it++)
    {
        ExMeasInfo item;

        // No conversion required for the seconds portion of the timestamps
        item.time_seconds = (*it).time_sec;

        // Convert the picoseconds fields to fractional seconds
        PICO pico = (PICO)(*it).time_pico_upper; // Don't shift until after we assign to avoid data loss
        pico <<= PICO_SHIFT;
        pico |= (*it).time_pico_lower;
        item.time_fraction = pico * PICOSECOND_EXPONENT;

        uint32_t payload_offset = 0;
        item.context_indicators.raw = (*it).context_indicator_field;

        if(item.context_indicators.ext_meas_info.trigger)
        {
            // No conversion required for the seconds portion of the timestamps
            item.trigger_seconds = (*it).context_fields[payload_offset];
            payload_offset++;

            // Convert the picoseconds fields to fractional seconds
            PICO trig_pico = (PICO)(*it).context_fields[payload_offset];
            payload_offset++;
            trig_pico <<= PICO_SHIFT;
            trig_pico |= (*it).context_fields[payload_offset];
            payload_offset++;
            item.trigger_fraction = trig_pico * PICOSECOND_EXPONENT;
        }

        if(item.context_indicators.ext_meas_info.ref_junction)
        {
            uint32_t ref_junction_word = (*it).context_fields[payload_offset];
            payload_offset++;

            typedef union
            {
                uint32_t u;
                float f;
            } convert;
            convert c;
            c.u = ref_junction_word;
            item.ref_junction = c.f;
        }

        retval.push_back(item);
    }
    return retval;
}

std::vector<DioData> convert_dio_data(std::vector<meas_fixed_point_dataset_t> data_vec, size_t samples_array_len)
{
    // Convert the data from the format used to store in HDF5 files to something that's a little easier to work with
    // This format is similar to the way data is returned from the Measurement.Read methods in the VTEXDigitizer/VTEXDsa drivers
    std::vector<DioData> retval;
    for(auto it = data_vec.begin(); it != data_vec.end(); it++)
    {
        DioData rec;

        // No conversion required for the seconds portion of the timestamp
        rec.time_seconds = (*it).time_sec;

        // Convert the picoseconds fields to fractional seconds
        PICO pico = (PICO)(*it).time_pico_upper; // Don't shift until after we assign to avoid data loss
        pico <<= PICO_SHIFT;
        pico |= (*it).time_pico_lower;
        rec.time_fraction = pico * PICOSECOND_EXPONENT;

        // Add the samples to a vector
        std::vector<int> sample_data;
        sample_data.assign(static_cast<uint32_t*>((*it).samples), static_cast<uint32_t*>((*it).samples) + samples_array_len);
        rec.samples = sample_data;

        retval.push_back(rec);
    }
    return retval;
}
