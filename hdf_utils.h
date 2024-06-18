#pragma once

#include <algorithm>
#include <ctime>
#include <iostream>
#include <map>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "H5Cpp.h"

// names for various groups, data sets, fields
const std::string root_group_name = "/";
const std::string ex_context_name = "EX_MEAS_INFO";
const std::string meas_float_name = "IF_MEAS_FLOAT32";
const std::string meas_fixed_point_name = "IF_MEAS_INT32";

const std::string time_seconds_field = "TimeSeconds";
const std::string time_pico_upper_field = "TimePicosecondsUpper";
const std::string time_pico_lower_field = "TimePicosecondsLower";
const std::string samples_field = "Samples";
const std::string context_indicator_field = "ContextIndicatorField";
const std::string context_words_field = "ContextFields";

// types and definitions for reading & processing data
#define MAX_SAMPLES 65536
#define MAX_CONTEXT_WORDS 128

typedef uint64_t PICO;
#define PICOSECOND_EXPONENT (1e-12)
#define PICO_SHIFT 32

#define MEAS_INT32_SHIFT 24

typedef struct meas_float_dataset_t {
    uint32_t time_sec;
    uint32_t time_pico_upper;
    uint32_t time_pico_lower;
    float samples[MAX_SAMPLES];
} meas_float_dataset_t;

typedef struct meas_fixed_point_dataset_t {
    uint32_t time_sec;
    uint32_t time_pico_upper;
    uint32_t time_pico_lower;
    uint32_t samples[MAX_SAMPLES];
} meas_fixed_point_dataset_t;

class RecordData {
public:
    int time_seconds;
    double time_fraction;
    std::vector<double> samples;
};

class RecordDataFloat {
public:
    int time_seconds;
    double time_fraction;
    std::vector<float> samples;
};

class DioData {
public:
    int time_seconds;
    double time_fraction;
    std::vector<int> samples;
};

typedef struct ex_meas_info_dataset_t {
    uint32_t time_sec;
    uint32_t time_pico_upper;
    uint32_t time_pico_lower;
    uint32_t context_indicator_field;
    uint32_t context_fields[MAX_CONTEXT_WORDS];
} ex_meas_info_dataset_t;

struct vrt_ext_context_indicators_meas_info_t {
    uint32_t unused : 24; // 0 - 23
    uint32_t ref_junction : 1; // 24
    uint32_t trigger : 1; // 25
    uint32_t comment : 1; // 26
    uint32_t eu : 1; // 27
    uint32_t weighting : 1; // 28
    uint32_t range : 1; // 29
    uint32_t span : 1; // 30
    uint32_t context_field_change_indicator : 1; // 31
};

typedef union {
    vrt_ext_context_indicators_meas_info_t ext_meas_info;
    uint32_t raw;
} vrt_context_presence_indicators_t;

class ExMeasInfo {
public:
    int time_seconds;
    double time_fraction;
    vrt_context_presence_indicators_t context_indicators;
    int trigger_seconds;
    double trigger_fraction;
    double ref_junction;

    ExMeasInfo() : time_seconds(0), time_fraction(0), trigger_seconds(0), trigger_fraction(0), ref_junction(0) {};
};
