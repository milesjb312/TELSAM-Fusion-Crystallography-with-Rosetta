#include <iostream>
using std::cout, std::endl, std::cerr, std::cin, std::getline;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <fstream>
using std::ifstream, std::ofstream;
#include <random>
using std::mt19937;
#include <stdio.h>

//To Do: remove neightboring low probability codons
//To Do: check for start codon followed by a HIS tag

map <char, vector<int>> ecoliProbMap = {
    {'F', {57, 100}},
    {'L', {15, 27, 39, 49, 54, 100}},
    {'I', {58, 93, 100}},
    {'M', {100}},
    {'V', {25, 43, 60, 100}},
    {'S', {11, 22, 37, 53, 67, 100}},
    {'P', {17, 30, 44, 100}},
    {'T', {16, 63, 76, 100}},
    {'A', {11, 42, 63, 100}},
    {'Y', {53, 100}},
    {'H', {55, 100}},
    {'Q', {30, 100}},
    {'N', {47, 100}},
    {'K', {73, 100}},
    {'D', {65, 100}},
    {'E', {70, 100}},
    {'C', {42, 100}},
    {'W', {100}},
    {'R', {36, 80, 87, 94, 96, 100}},
    {'G', {29, 75, 88, 100}},
    {'*', {0, 64, 100}},
};

map <char, vector<string>> codonMap = {
    {'F', {"TTT", "TTC"}},
    {'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}},
    {'I', {"ATT", "ATC", "ATA"}},
    {'M', {"ATG"}},
    {'V', {"GTT", "GTC", "GTA", "GTG"}},
    {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},
    {'P', {"CCT", "CCC", "CCA", "CCG"}},
    {'T', {"ACT", "ACC", "ACA", "ACG"}},
    {'A', {"GCT", "GCC", "GCA", "GCG"}},
    {'Y', {"TAT", "TAC"}},
    {'H', {"CAT", "CAC"}},
    {'Q', {"CAA", "CAG"}},
    {'N', {"AAT", "AAC"}},
    {'K', {"AAA", "AAG"}},
    {'D', {"GAT", "GAC"}},
    {'E', {"GAA", "GAG"}},
    {'C', {"TGT", "TGC"}},
    {'W', {"TGG"}},
    {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},
    {'G', {"GGT", "GGC", "GGA", "GGG"}},
    {'*', {"TAG", "TAA", "TGA"}},
};

map<string, string> improbableCodons = {
    {"TTA", "CTG"},{"TTG", "CTG"},
    {"CTT", "CTG"},{"CTC", "CTA"},
    {"CTA", "CTG"},{"ATA", "ATT"},
    {"TCT", "AGC"},{"TCC", "AGC"},
    {"TCA", "AGC"},{"CCC", "CCG"},
    {"CCA", "CCG"},{"ACA", "ACC"},
    {"GCT", "GCG"},{"CGA", "CGC"},
    {"CGG", "CGC"},{"AGT", "AGC"},
    {"AGA", "AGC"},{"AGG", "AGC"},
    {"GGA", "GGC"},{"GGG", "GGC"}
};


string readFasta(ifstream& inputstream){
    cout << "inside readFasta function" << endl;
    if (!inputstream){
        throw std::invalid_argument("could not open file");
    }
    string line;
    while (getline(inputstream, line)) {
        if (!line.empty() && line[0] != '>') {
            break;
        }
    }
    return line;
}


//bool requestHisTag(){
//    string response;
//    cout << "Do you want to add a His-tag? [Y/N]" << endl;
//    getline(cin, response);
//    if(response == "Y" || response == "y"){
//        return true;
//    }
//    return false;
//}

void prepSequence(string& sequence){
    for (char& c : sequence){
        c = toupper(c);
        if (codonMap.find(c) == codonMap.end()){
            cerr << "found invalid sequence character: " << c << endl;
            throw std::invalid_argument("invalid amino acid");
        }
    }
    
    //if(requestHisTag()){
    //    cout << "Adding a 10x His-tag" << endl;
    //    sequence = "HHHHHHHHHH" + sequence;
    //}
    if(sequence.front() != 'M'){
        cout << "Adding a start codon" << endl;
        sequence = "M" + sequence;
    }
    if(sequence.back() != '*'){
        cout << "Adding stop codon" << endl;
        sequence = sequence + "*";
    }
}

void checkPairs(string& sequence){
    int i = 1;
    string prev = sequence.substr(0,3);
    string current;
    string newCodon;
    while (i < sequence.size()/3){
        current = sequence.substr(i*3, 3);
        //cout << "Index: " << i*3 << " Previous: " << prev << " Current: " << current << endl;
        if (improbableCodons.find(prev) != improbableCodons.end() && improbableCodons.find(current) != improbableCodons.end()){
            newCodon = improbableCodons[current];
            sequence.replace(i*3, 3, newCodon);
            cout << "Replacing " << current << " with " << newCodon << " at index " << i*3 << endl;
            current = newCodon;
        }
        prev = current;
        i++;
    }
    return;
}

int seed = 100;
mt19937 gen(seed);
std::uniform_int_distribution<> dist(1,100);

string selectCodon(char AA){
    int index = 0;
    int rand_selector = dist(gen); 
    vector<int> list = ecoliProbMap[AA];
    while(index < list.size() && list.at(index) < rand_selector){
        index++;
    }
    return codonMap[AA].at(index);
}

string reverseTranslate(string AAseq){
    string DNAseq = "";
    for (char c : AAseq){
        DNAseq = DNAseq + selectCodon(c);
    }
    return DNAseq;
}

int main(int argc, char* argv[]){

    string inputSequence;
    ifstream infile(argv[1]);
    inputSequence=readFasta(infile);
    infile.close();

    // Now we need to check the sequence and prepare it for reverse translation
    if(inputSequence==""){return -1;}

    cout << "This is the input sequence being used:" <<endl;
    cout << inputSequence << endl;

    prepSequence(inputSequence);

    cout << inputSequence << endl;

    // Now we reverse translate
    string outputSequence = reverseTranslate(inputSequence);

    checkPairs(outputSequence);

    //print sequence to output folder
    string outfileName;
    outfileName = argv[2];
    ofstream outfile(outfileName);
    string header;
    header = argv[3];

    if (!outfile.is_open())
    {
        cerr <<"Could not open could not open the file but I will print it to the screen for you :)" << endl;
        cout << outputSequence << endl;
        return 2;
    }

    outfile << "> " << header << endl;
    outfile << outputSequence << endl;

    outfile.close();

    cout << "DONE! \nClosing the program" << endl;
    return 0;
    
}