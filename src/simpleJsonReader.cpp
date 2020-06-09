#include "simpleJsonReader.hpp"
// Loads debug_settings structure from the specified XML file
void SimpleJsonReader::load(const std::string &filename)
{
    // Load the json file into the property tree. If reading fails, exception raised
    read_json(filename,pt);
}

std::string SimpleJsonReader::extractValue(const std::string& key){
    ptree& tree = pt.get_child(key); 
    std::string type = tree.get<std::string>("type");
    std::string tmp =  tree.get<std::string>("value");
    return tmp;
}

std::vector<std::string> SimpleJsonReader::extractArray1D(const std::string& key){
    ptree& tree = pt.get_child(key); 
    std::vector<std::string> tmp;
    BOOST_FOREACH(ptree::value_type &v,tree.get_child("value")){    
        tmp.push_back(v.second.data());
    } 
    return tmp;
}
std::vector<std::vector<std::string>> SimpleJsonReader::extractArray2D(const std::string& key){
    ptree& tree = pt.get_child(key); 
    std::string type = tree.get<std::string>("type");
    std::vector<std::vector<std::string>> tmp;
        BOOST_FOREACH(ptree::value_type &layer1,tree.get_child("value")){    
            std::vector<std::string> array1D;
            BOOST_FOREACH(ptree::value_type &layer2,layer1.second){
              array1D.push_back(layer2.second.data());   
            } 
            tmp.push_back(array1D);
        }
    return tmp;
}