#ifndef SIMPLEJSONREADER_HPP
#define SIMPLEJSONREADER_HPP
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <iostream>
using boost::property_tree::ptree;
struct SimpleJsonReader
{
    ptree pt;
    std::string extractValue(const std::string& key);
    std::vector<std::string> extractArray1D(const std::string& key);
    std::vector<std::vector<std::string>> extractArray2D(const std::string& key);
    void load(const std::string &filename);
};
#endif