#ifndef FILEREADINGMETHODS_HXX
#define FILEREADINGMETHODS_HXX

#include<string>
#include<vector>

namespace disjointPaths {


template<class T=char>
 std::vector<std::string> split(
                std::string inputString, T delim) {
        size_t occurence = 0;
        size_t newOccurence = 0;
        std::vector<std::string> strings;
        while (newOccurence < inputString.size()) {
                newOccurence = std::min(inputString.find_first_of(delim, occurence),
                                inputString.size());

                std::string newString(inputString, occurence, newOccurence - occurence);
                strings.push_back(newString);
                newOccurence = newOccurence + 1;
                occurence = newOccurence;
        }

        return strings;
}

}


#endif // FILEREADINGMETHODS_HXX
