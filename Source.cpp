#include <iostream>    
#include <vector>               
#include <utility>                  
#include <algorithm>     
#include <fstream>   
#include <boost/range/irange.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/regex.hpp>
#include <boost/config.hpp>
#include<ctime>
#include<chrono>
using namespace std::chrono;
#include <math.h>

using namespace boost;
using namespace std;

int const n_attr = 4; // number of attributes
int const a_max = -2; // value for non existing attributes
struct Vertex
{
    int Attribute[n_attr]; // matrix of each value's attribute Attribute[0] = a1, Attribute[1] = a2, etc etc... 
};
// Define the structure of graph
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Vertex> Graph;

Graph read_graph(string edges_filename, string att_filename) 
{
    /// <summary>
    /// Function to read raw graph data
    /// </summary>
    /// <param name="edges_filename"> text files containing edges
    /// <param name="att_filename"> text file containing vertex attributes
    /// <returns> filled graph

    Graph g; // Final graph
    string line_edges, line_att;
    int nb_line_att = 0;
    ifstream file_edges, file_attributes, file_label;
    
    file_attributes.open(att_filename);
    if (file_attributes.is_open()) // Check if the file is open
    {
        cout << "File " + att_filename + " opened successfully" << endl;
    }
    vector <string> Att;
    while (file_attributes.good()) // We start to read the file
    {                              //Remark: using the loop while(file.good()) will make a read until a line past the last non empty line
        getline(file_attributes, line_att); // Read line by line.
        split(Att, line_att, is_any_of("\t")); // Splits the line into a vector containing all strings of a text line based on their separation (in this case "\t")
        if (Att.size() > 1) // Condition for reading non empty lines. This may need further modification if we have a single attribute per line
        {
            add_vertex(g); // Create a vertex
            // In the following loop we write the attributes
            for (int i = 0; i < Att.size(); i++) // We use .size()-1 as the split function creates a vector of size n + 1
                                                    // (this can be set up according to the nature of the file we read)
            {
                int att;
                // Transform strings into integers
                istringstream(string(Att[i])) >> att;
                g[nb_line_att].Attribute[i] = att;
            }
            nb_line_att++; // nb_line_att refers to the index of a vertex in the graph
        }
    }
    file_attributes.close(); // Close file
    
    file_edges.open(edges_filename);
    if (file_edges.is_open())
    {
        cout << "File " + edges_filename + " opened successfully" << endl;
    }
    vector <string> Edges;
    while (file_edges.good())
    {
        getline(file_edges, line_edges); // Read line
        split(Edges, line_edges, is_any_of("	"));
        if (Edges.size() == 2) { // Check if the line is not empty (otherwise the size would be 1)
            int ledge, redge;
            // Transform strings into integers
            istringstream(string(Edges[0])) >> ledge;
            istringstream(string(Edges[1])) >> redge;
            //for (auto i: vertices(g))

            add_edge(ledge, redge, g); // Add edges
        }
    }
    file_edges.close();

    return g;
}

vector<vector<int>> CombineAttributesSize(vector<int>  A, vector<int> D, int k)
{
    /// <summary>
    /// Create combination attributes of size k
    /// </summary>
    /// <param name="A"> set of attributes
    /// <param name="D"> set of attribute values
    /// <param name="k"> combination size
    /// <returns> combination attributes of size k

    vector<vector<int>> attcombint; // final vector
    for (int i = 0; i < A.size() - k + 1; i++) // Create combination of size k.
    {
        for (int j = 0; j < D.size(); j++)
        {
            vector< vector<int>> attcomb;
            attcomb.push_back({ A[i],D[j] });  // Add size 1 combination 
            int m = 1;
            while (m < k) // If k = 1, doesn't enter in the loop
            {
                vector<vector<int>> newcomb;

                for (vector<int> elem : attcomb)
                {
                    for (int ind = elem[elem.size() - 2] + 1; ind < A.size(); ind++) // Add 1 to move to the next attribute
                    
                    {
                        for (int j = 0; j < D.size(); j++) // Extend the combination with next size 1 combination
                        {
                            vector<int> inter = elem;
                            inter.push_back(A[ind]);
                            inter.push_back(D[j]);
                            newcomb.push_back(inter);

                        }
                    }
                }
                attcomb = newcomb; // Save the combination in construction
                m++;
            }
            attcombint.insert(attcombint.end(), attcomb.begin(), attcomb.end()); // Extraction of each combination
        }
    }
    return attcombint;
}

vector<vector<int>> CombineAttributes(vector<int>  A, vector<int> D) 
{
    /// <summary>
    /// Construct all combinations of attributes
    /// </summary>
    /// <param name="A"> set of attributes
    /// <param name="D"> set of attribute values
    /// <returns> combination attributes

    vector< vector<int> > attcomb; // Final vector
    for (int k = 1; k <= A.size(); k++)
    {
        vector < vector<int>> temp = CombineAttributesSize(A, D, k);
        attcomb.insert(attcomb.end(), temp.begin(), temp.end());
    }
    return attcomb;
}

vector<vector<int>> CreatedConnectedComponents(vector<vector<int>> AD, int minvol, int maxvol)
{
    /// <summary>
    /// Create all pattern possibilities respecting volume constraint, for one graph.
    /// </summary>
    /// <param name="AD"> Attribute combination set
    /// <param name="minvol"> Minimum size contraint
    /// <param name="maxvol"> Maximum size constraint
    /// <returns> Vector pattern possibilities

    vector<vector<int> > connectcomp; // Final vector
    int m = minvol;
    for (int i = 0; i <= maxvol - minvol; i++) // Create combination of minvol-maxvol.
    {
        for (int j = 0; j < AD.size(); j++)
        {
            vector< vector<int>> temp;
            temp.push_back(AD[j]);  // Add first vertex

            for (int n = 1; n < m; n++) // If minvol = 1, doesn't enter in the loop, m travels from minvol to maxvol
            {
                vector<vector<int>> newcomb;
                for (vector<int> elem : temp)
                {
                    for (vector<int> z : AD) // Extend the combination with next vertex
                    {
                        vector<int> inter = elem;
                        inter.push_back(-3); // Separation between vertex
                        for (int k = 0; k < z.size(); k++)
                        {
                            inter.push_back(z[k]);
                        }
                        newcomb.push_back(inter);
                    }
                }
                temp = newcomb; // Save the vector in construction
            }
            connectcomp.insert(connectcomp.end(), temp.begin(), temp.end());
        }
        m++;
    }
    return connectcomp;
}

vector<vector<int>> CombineTimeSize(vector<int> T, int k)
{
    /// <summary>
    /// Create combination time starting by t1 of size k
    /// </summary>
    /// <param name="T"> Set of time
    /// <param name="k"> Combination size
    /// <returns> Combination time starting by t1

    vector<vector<int>> timecombint, timecomb; // Final vector
    timecombint.push_back({ T[0] });  // Add size 1 combination 
    int m = 1;
    while (m < k) // If k = 1, doesn't enter in the loop
    {
        vector<vector<int>> newcomb;
        for (vector<int> elem : timecombint)
        {
            for (int ind = elem[elem.size() - 1] + 1; ind < T.size(); ind++) // Need T start by 0 to work well
            {
                vector<int> inter = elem;
                inter.push_back(T[ind]);
                newcomb.push_back(inter);

            }
        }
        timecombint = newcomb;
        m++;
    }
    timecomb.insert(timecomb.end(), timecombint.begin(), timecombint.end()); // extraction of each combination
    return timecomb;
}

vector<vector<int>> CombineTime(vector<int>  T) // show all combinations of times starting by 1
{
    /// <summary>
    /// Construct all combination time
    /// </summary>
    /// <param name="T"> Set of time
    /// <returns></returns>

    vector<vector<int>> timecomb; // Final vector
    for (int k = 1; k <= T.size(); k++)
    {
        vector<vector<int>> temp = CombineTimeSize(T, k);
        timecomb.insert(timecomb.end(), temp.begin(), temp.end());
    }
    return timecomb;
}
// End combinatories functions

// Combinatories Graph Function
vector<vector<int>> Occurence2(Graph g, vector<int> mot) // return the vector of vertex who's respect the motif
{
    /// <summary>
    /// Search in a graph all occurrences that follow a pattern
    /// </summary>
    /// <param name="g"> Graph
    /// <param name="mot"> Pattern to find
    /// <returns> vector of vertices that respect the pattern

    boost::graph_traits <Graph>::vertex_iterator i, end; // Create iteration on vertices
    vector<vector<int>> occ; // Final vector
    vector <int> vertexvis;
    for (tie(i, end) = vertices(g); i != end; ++i) // Explore all vertices of graph
    {
        vector<vector<int>> comb;
        int m = 0, val = *i;
        bool test = true;
        while (m < mot.size() and mot[m] != -3) // Condition of end of pattern or end of vertex condition
        {
            if (g[*i].Attribute[mot[m]] != mot[m + 1]) // Check if vertex attribute in the graph respect the attributes in the pattern 
            {
                test = false;
                break;
            }
            m += 2;
        }
        if (test)
        {
            comb.push_back({ val });
            if (m < mot.size())
            {
                if (mot[m] == -3) // Check separation vertex
                {
                    m++;
                }
            }
            while (m < mot.size()) // Check the end of pattern
            {
                if (mot[m] == -3)
                {
                    m++;
                }
                int n = 0;

                vector<vector<int>> newcomb;

                for (vector<int> temp : comb)
                {

                    auto neighbours = adjacent_vertices(temp[temp.size() - 1], g); // Function to find neigbours of vertex

                    for (size_t vd : make_iterator_range(neighbours)) // explore all neighbours of vertices
                    {
                        bool testvertex = true;

                        for (int elem : vertexvis) // check if vertex is already visited in the graph 
                        {
                            if (elem == vd)
                            {
                                testvertex = false;
                            }
                        }

                        if (testvertex)
                        {

                            test = true;
                            for (int elem : temp) // check if vertex is already in the pattern
                            {
                                if (elem == vd)
                                {
                                    test = false;
                                    break;
                                }
                            }

                            if (test) // Vertex not in the pattern
                            {
                                bool test2 = true;
                                n = m;
                                while (n < mot.size() and mot[n] != -3) // Check separation of vertex and end of pattern
                                {
                                    if (g[vd].Attribute[mot[n]] != mot[n + 1]) // Check if attribute's vertex respect the pattern
                                    {
                                        test2 = false;
                                        break;
                                    }
                                    n += 2;
                                }

                                if (test2) // Keep all sequence element who's respect the condition and we add it to the solution
                                {
                                    vector<int> inter = temp;
                                    inter.push_back(vd);
                                    newcomb.push_back(inter);

                                }

                            }
                        }
                    }
                }
                comb = newcomb; // Update the vector of occurrence in construction
                n = m;
                while (n < mot.size() and mot[n] != -3) // Run throught the next vertex in the pattern
                {
                    n += 2;
                }
                m = n;
            }
        }
        if (comb.size() > 0)
            occ.insert(occ.end(), comb[0]);
        vertexvis.push_back(val); // Add the vertex visited 
    }
    if (occ.size() < 1) // Check the condition minvol and if a graph contains the motif
    {
        occ.push_back({ {-1} });
    }
    return occ;
}

vector<vector<vector<int>>> OccurenceTimeCombination2(vector<vector<vector<int>>>  setocc) 
{
    /// <summary>
    /// Addition occurences of graph with other graphs respecting pattern
    /// </summary>
    /// <param name="setocc"> Set of occurences
    /// <returns> Union of occurence separate by -5 to difference the time graph
    vector<vector<vector<int>>> combsol; // Final vector
    vector<vector<int>> combsolinter = setocc[0];
    for (int i = 1; i < setocc.size(); i++)
    {
        combsolinter.push_back({ -5 }); // Separation of sequence vertices
        combsolinter.insert(combsolinter.end(), setocc[i].begin(), setocc[i].end());
    }
    combsol.push_back(combsolinter);
    return combsol;
}

std::tuple<vector<vector<int>>, vector<vector<vector<int>>>, vector<int>, vector<vector<int>>, vector<int>> Extension
    (vector<vector<vector<int>>> combocc1, vector<vector<vector<int>>> combocc2,
    vector<vector<int>>combmot1, vector<vector<int>>combmot2, vector<int> count1, vector<int> count2,
    int mincom, int minsup)
{
    /// <summary>
    /// Extend the pattern with a new pattern if they have enought commun vertices (mincom) and minimum frequency, count the new frequency, 
    /// keep the extension of the occurence
    /// </summary>
    /// <param name="combocc1"> Set of occurence to check the extension
    /// <param name="combocc2"> Set of occurence to check extend with
    /// <param name="combmot1"> Set of pattern to extend
    /// <param name="combmot2"> Set of pattern extend with
    /// <param name="count1"> Set of frequency 
    /// <param name="count2"> Set of frequency
    /// <param name="mincom"> Minimum commum vertex constraint
    /// <param name="minsup"> Minimum frequency constraint
    /// <returns> Extended Pattern, non extended pattern, set of frequency, final frequency (for non extended ), new set of occurence

    vector<vector<vector<int>>> occextend; // result extension occurence
    vector<vector<int>> extendsolmot, notextendsolmot; // result extend motif and final motif
    vector<int> countextend, countfinal; // result extend frequency and final frequency
     for (int occ1 = 0; occ1 < combocc1.size(); occ1++) 
    {
        vector<vector<int>> extendsolmotinter;

        for (int occ2 = 0; occ2 < combocc2.size(); occ2++)
        {
            vector<vector<int>> tempocc;
            vector<int> sep = {-5};
            auto it = find(combocc2[occ2].begin(), combocc2[occ2].end(), sep); 

            int lastocc = 0, pos = 0, pos2 = it - combocc2[occ2].begin(); // return position of -5 if exist else last element of list
            for (vector<int>elem : combocc1[occ1]) // run throught the first set of occurence
            {
                if (elem[0] == -5)
                {
                    tempocc.push_back(elem);
                    lastocc = tempocc.size();
                    pos = pos2+1;
                    auto it = find(combocc2[occ2].begin()+pos, combocc2[occ2].end(), elem);
                    if (it == std::end(combocc2[occ2]))
                    {
                        pos2 = combocc2[occ2].size();
                    }                       
                    else
                    {
                        pos2 = it - combocc2[occ2].begin();
                    }
                    
                }
                else
                {
                    for (int elem2 = pos; elem2 < pos2; elem2++)
                    {
                        int minc = 0;
                        for (int vertex : elem) // We keep all values
                        {
                            for (int vertex2 : combocc2[occ2][elem2])
                            {
                                
                                if (vertex == vertex2)
                                {
                                    minc++;
                                }
                            }
                        }
                        if (minc >= mincom)
                        {
                            
                            if(tempocc.size() > 0)
                            {
                                auto occ = find(tempocc.begin()+lastocc, tempocc.end(), combocc2[occ2][elem2]);

                                if (occ == std::end(tempocc))
                                {
                                    tempocc.push_back(combocc2[occ2][elem2]);
                                }

                            } 
                            else
                            {
                                tempocc.push_back(combocc2[occ2][elem2]);
                            }                      
                            
                        }
                    }
                }      
            }
        
            if(tempocc.size()>=minsup)
            {
                vector<int> mot = combmot1[occ1];
                int count = tempocc.size();
                mot.push_back(-5);
                mot.insert(mot.end(), combmot2[occ2].begin(), combmot2[occ2].end());
                extendsolmotinter.push_back(mot); // add next motif to the pattern
                countextend.push_back(count); // update new frequency
                occextend.push_back(tempocc); // keep the vertices to extend
            }
        }

        extendsolmot.insert(extendsolmot.end(), extendsolmotinter.begin(), extendsolmotinter.end());

        if (extendsolmotinter.size() < 1)
        {
            notextendsolmot.push_back(combmot1[occ1]); // Final pattern
            countfinal.push_back(count1[occ1]); // Final frequency
        }
    }
    return  make_tuple(extendsolmot, occextend, countextend, notextendsolmot, countfinal);
}

vector<vector<vector<vector<int>>>> CreatePatterns(vector<Graph> gt, vector<int> A, vector <int> D, int minvol, int maxvol, int minsup,int maxsup, int mincom)
{ // Pensez à enlever maxsup lors du push de l'algo

    /// <summary>
    /// Search on temporal attributed graph patterns repecting constraints
    /// </summary>
    /// <param name="gt"> Vector of graphs
    /// <param name="A"> Vector of attributes
    /// <param name="D"> Vector of attributes values
    /// <param name="minvol"> minvol constraint
    /// <param name="maxvol"> maxvol constraint
    /// <param name="minsup"> minsup constraint
    /// <param name="mincom"> mincom constraint
    /// <returns> All patterns on graphs respecting min-maxvol, minsup and mincom constraints
    // Initialization of pattern
    int num_attr = A.size(), num_time = gt.size(), cpt = 1;
    vector<int> time; push_back(time, irange(0, num_time));
    vector<vector<int>> combatt = CombineAttributesSize(A, D, num_attr); 
    vector<vector<int>> compconnect = CreatedConnectedComponents(combatt, minvol, maxvol);
    vector<vector<int>> combtime = { {} };
    vector<vector<vector<vector<int>>>> Result; // Final vector

    vector<vector<vector<vector<int>>>> occurence; // Vector of occurences in each graph


    for (int i : time) // loop to search for patterns in graphs
    {
        vector<vector<vector<int>>> occinter;
        auto t_start = std::chrono::high_resolution_clock::now();

        for (vector<int> motif : compconnect) // run throught all motif between minvol, maxvol
        {
            occinter.push_back(Occurence2(gt[i], motif));
        }
        occurence.push_back(occinter);
    }
    for (vector<int> t1 : combtime) // Run throught all t1 combination
    {
        vector<int> t1comb = t1;
        int lastt1comb = t1comb[t1comb.size() - 1];
        vector<int> count, veccount;
        vector<vector<int>> seqmotif, seqmotiftime;
        vector<vector<vector<int>>> seqinter, seqminocc;
        //////// Initialisation

        if (t1comb.size() < 2) // separate case if we take a singleton of time or sets of times to save calculation

        {
            for (int mot = 0; mot < compconnect.size(); mot++)
            {
                if (occurence[0][mot].size() > 1)
                {
                    int min = occurence.size();
                    if (maxsup >min and min >= minsup) // Check minsup constraint
                    {

                        seqmotif.push_back(compconnect[mot]);
                        count.push_back(min);
                        seqinter.push_back(occurence[0][mot]);
                    }
                }
            }
        }

        else
        {

            for (int mot = 0; mot < compconnect.size(); mot++) // mot is the line index of a pattern
            {
                int min = 0;

                for (int i : t1comb)
                {
                    if (occurence[i][mot][0].size() > 1)
                    {
                        min += occurence[i][mot].size();
                    }
                }

                if (maxsup > min and min >= minsup) // Check minsup constraint for multiple graph
                {
                    vector<vector<vector<int>>> minocc;
                    for (int i : t1comb)
                    {
                        minocc.push_back(occurence[i][mot]);

                    }
                    minocc = OccurenceTimeCombination2(minocc); // Create additioned graph                                              

                    seqmotif.push_back(compconnect[mot]); // Inter pattern 
                    count.push_back(min); // Inter frequency
                    seqinter.insert(seqinter.end(), minocc.begin(), minocc.end()); // Inter occurrence

                }
            }
        }
        //// Intermediate values
        veccount = count;
        seqminocc = seqinter;
        seqmotiftime = seqmotif;

        for (int t = 1; t < time.size() - lastt1comb; t++) // Run throught all graph depend of the last elem of the initialisation
        {

            vector<vector<int>> seqmotif, seqmotif2;
            vector<vector<vector<int>>> seqinter;
            vector<int> count;
            vector<int> countfinal;
            auto t_start2 = std::chrono::high_resolution_clock::now();
            if (t1comb.size() < 2) // Separate case if we take a singleton of time or sets of times to save calculation
            {
                for (int mot = 0; mot < compconnect.size(); mot++) // mot is the line index of a pattern
                {
                    if (occurence[0][mot][0].size() > 1)
                    {
                        int min = occurence[0][mot].size();
                        if ( min  >= minsup) // Check minsup constraint 
                        {

                            seqmotif.push_back(compconnect[mot]);
                            count.push_back(min);
                            seqinter.push_back(occurence[0][mot]);
                        }
                    }
                }
            }

            else
            {

                for (int mot = 0; mot < compconnect.size(); mot++) // mot is the line index of a pattern
                {
                    int min = 0;

                    for (int i : t1comb)
                    {
                        if (occurence[i][mot][0].size() > 1)
                        {
                            min += occurence[i][mot].size();
                        }
                    }

                    if (maxsup>min and min >= minsup) // Check minsup constraint
                    {
                        vector<vector<vector<int>>> minocc;
                        for (int i : t1comb)
                        {
                            minocc.push_back(occurence[i][mot]);

                        }
                        minocc = OccurenceTimeCombination2(minocc); // Create additioned graph

                        seqmotif.push_back(compconnect[mot]);
                        count.push_back(min);
                        seqinter.insert(seqinter.end(), minocc.begin(), minocc.end());

                    }
                }
            }
            /////// Extension Part

            auto extinter = Extension(seqminocc,
                seqinter, seqmotiftime, seqmotif, veccount, count, mincom, minsup);
            seqmotiftime = get<0>(extinter); seqminocc = get<1>(extinter); veccount = get<2>(extinter);
            seqmotif2 = get<3>(extinter); countfinal = get<4>(extinter);

            if (seqmotif2.size() > 0) // Result for none extend pattern
            {

                for (int i = 0; i < seqmotif2.size(); i++)
                {
                    vector<vector<vector<int>>> temp;
                    temp.push_back({ t1comb,{time[t]} });
                    temp.push_back({ {countfinal[i]} });
                    temp.push_back({ seqmotif2[i] });
                    Result.push_back(temp);
                }

            }

            if ((seqmotiftime.size() < 1) and (t + 1 < time.size() - lastt1comb)) // If there are no more patterns to extend
            {// Start a new pattern on the next time
                t++;
                for (int motif = 0; motif < compconnect.size(); motif++) 
                {

                    int min = 0;


                    for (int i : t1comb)
                    {
                        if (occurence[i + t][motif][0].size() > 1)
                        {
                            min += occurence[i + t][motif].size();
                        }
                    }


                    if ( maxsup > min and min  >= minsup)
                    {
                        vector<vector<vector<int>>> minocc;
                        for (int i : t1comb)
                        {
                            minocc.push_back(occurence[i + t][motif]);
                        }
                        minocc = OccurenceTimeCombination2(minocc);

                        seqmotif.push_back(compconnect[motif]);
                        count.push_back(min);
                        seqinter.insert(seqinter.end(), minocc.begin(), minocc.end());

                    }

                }
                veccount = count;
                seqminocc = seqinter;
                seqmotiftime = seqmotif;

            }

        }
        if (seqmotiftime.size() > 0)
        {
            for (int i = 0; i < seqmotiftime.size(); i++)
            {
                vector<vector<vector<int>>> temp;
                temp.push_back({ t1comb,{time[time.size() - 1] + 1} });
                temp.push_back({ {veccount[i]} });
                temp.push_back({ seqmotiftime[i] });
                Result.push_back(temp);
            }
        }

    }
    return Result;
}


int main()
{
    vector<Graph> gt;
    vector<vector<vector<vector<int>>>> result;
    
     //Import CovidChina
     for (int i = 0; i < 2; i++)
    {
        Graph g;
        stringstream stmp2;
        stringstream stmp;
        stmp2 << "./data/covid_china/data/attributes/T" << std::to_string(i) << "attributes";
        stmp << "./data/covid_china/data/attributes/T" << std::to_string(i) << "edges";
        g = read_graph(stmp.str(), stmp2.str(), {0,1,2});
        gt.push_back(g);
    }
    //int minsup =260*.15 , maxsup = 100000;
    int minsup = 10, maxsup = 1000;
    CreatePatterns(gt,{0,1,2},{1,0,2,3},3,3, minsup, maxsup, 3);
    return 0;
}
