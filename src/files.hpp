
#ifndef FILES_HPP_
#define FILES_HPP_

/** Create directory
 *
 * All subdirectories are created as necessary.
 * */
void create_directory(const std::string& dir);

/** Remove file or directory
 *
 * This is equivalent to 'rm -rf fname'
 * */
void remove_file(const std::string& fname);

/** Compress file iname to oname.zip using zip
 *
 * The file is moved to the archive, i.e. iname does not exsit after this
 * function has returned.
 * */
void compress_file(const std::string& iname, const std::string& oname);

#endif//FILES_HPP_
