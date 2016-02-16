/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __MISC_H__
#define __MISC_H__

#include <istream>
#include <vector>
#include <string>

void read_unitigs_lengths(std::istream& is, std::vector<int>& lengths);
void read_unitigs_sequences(std::istream& is, std::vector<int>& lengths);
void read_unitigs_sequences(std::istream& is, std::vector<int>& lengths,
                            std::vector<std::string>& sequences);


#endif /* __MISC_H__ */
