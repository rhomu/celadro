#include "header.hpp"
#include "model.hpp"
#include "simpleJsonReader.hpp"
#include <exception>
#include <boost/lexical_cast.hpp>
using namespace std;
inline bool isExisted(const string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}
template<class T,size_t D>
inline vec<T,D> convert2Vec(const std::vector<std::string>& strArr){
  vec<T,D> tmp;
  for (unsigned i=0;i<D;i++){
    tmp[i] = boost::lexical_cast<T>(strArr[i]);
  }
  return tmp;
}
template<class T,size_t D>
inline vector<vec<T,D>> convert2VectorOfVec(const vector<vector<string>>& strArr){
  vector<vec<T,D>> tmp;
  for(auto it1 = strArr.begin(); it1 != strArr.end(); it1++){
    vec<T,D> vec0;
    for (unsigned i = 0; i < D;i++){
        vec0[i] = boost::lexical_cast<T>((*it1)[i]);  
    }
    tmp.push_back(vec0);
  }
  return tmp;
}
template<class T>
inline vector<T> convert2VectorOfPtype(const vector<string>& strArr){
  vector<T> tmp;
  for(auto it1 = strArr.begin(); it1 != strArr.end(); it1++){
    T t = boost::lexical_cast<T>(*it1);  
    tmp.push_back(t);
  }
  return tmp;
}
void Model::ReadConfig(const std::string& input_dir_file){
    if (!isExisted(input_dir_file)){
      throw error_msg("Error when opening initial configuration .json file");
    }
    SimpleJsonReader sjr; 
    try{
       sjr.load(init_config_file);
    }
    catch(std::exception e){
       std::cerr<<e.what()<<std::endl;
       throw error_msg("Error when loading.json file");
    }
    unsigned nphases0 = std::stoi(sjr.extractValue("nphases"));
    unsigned BC0 = std::stoi(sjr.extractValue("BC"));
    if ((nphases != nphases0) || (BC !=BC0)){
      throw error_msg("Mismatch between init_config and runcard(Boundary condition and nphases mismatch)");
    }
    //init_cell_shape = std::stoi(sjr.extractValue("init_cell_shape"));
    //init_aspect_ratio = std::stoi(sjr.extractValue("init_aspect_ratio"));
    //vector<string> tmp = sjr.extractArray1D("Size");
    vec<unsigned,2> Size0 =  convert2Vec<unsigned,2>(sjr.extractArray1D("Size"));
    if ((Size0[0]!=Size[0]) || (Size0[1]!=Size[1])){
      throw error_msg("Mismatch between init_config and runcard(Size mismatch)");
    }
    std::vector<vector<string>> tmp = sjr.extractArray2D("init_centers");
    init_centers = convert2VectorOfVec<unsigned,2>(tmp);
    theta_pol = convert2VectorOfPtype<double>(sjr.extractArray1D("theta_pol"));
    theta_nem = convert2VectorOfPtype<double>(sjr.extractArray1D("theta_nem"));
    alpha = convert2VectorOfPtype<double>(sjr.extractArray1D("alpha"));
}
