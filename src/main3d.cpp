#include <iostream>
#include <fstream>
#include <random>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <cstdlib>
#include<vector>
#include<algorithm>
#include<numeric>
#include<sstream>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <prettyprint.hpp>

#define DEBUG (false)
namespace po=boost::program_options;
  

int L;
int N;

std::vector<int> sites;
std::vector<int> parent;

void init(int size){
    L=size;
    N = L * L * L;
    parent.resize(N);
    sites.resize(N);
    for(int i = 0; i< N; ++i ){
        parent[i] = i ;
    }
}
    
int find(int i){
    if (i != parent[i]){
        parent[i] = find(parent[i]); // tree is flattened; https://algo-logic.info/union-find-tree/
    }
    return parent[i] ;
}

void unite(int i, int j){  //can be optimized more if necessary 
    i = find(i);
    j = find(j);
    parent[j] = i ;  
}

void connect(int i, int j){
    if(sites[i] == 0){
        return ;
    }
    if(sites[j] == 0){
        return ;
    }
    unite(i,j);
}

int pos2index(int ix, int iy, int iz){
    ix =(ix + L) % L;
    iy =(iy + L) % L;
    iz =(iz + L) % L;
    return ix *L*L + iy *L + iz ;
}

double crossing_probability_z(void){
    for(int ix1 = 0; ix1 < L ; ix1++){
    for(int iy1 = 0; iy1 < L ; iy1++){
        int i = pos2index(ix1, iy1, 0);
        int ci = find(i);
        for(int ix2 = 0 ; ix2 < L; ++ix2){
        for(int iy2 = 0 ; iy2 < L; ++iy2){
            int j = (pos2index(ix2, iy2, L-1));
            int cj = find(j);
            if(ci == cj){
                return 1.0;
            }
        }
        }
    }
    }
    return 0.0;
}

double crossing_probability_y(void){
    for(int ix1 = 0; ix1 < L ; ix1++){
    for(int iy1 = 0; iy1 < L ; iy1++){
        int i = pos2index(ix1, 0, iy1);
        int ci = find(i);
        for(int ix2 = 0 ; ix2 < L; ++ix2){
        for(int iy2 = 0 ; iy2 < L; ++iy2){
            int j = (pos2index(ix2, L-1, iy2));
            int cj = find(j);
            if(ci == cj){
                return 1.0;
            }
        }
        }
    }
    }
    return 0.0;
}

double crossing_probability_x(void){
    for(int ix1 = 0; ix1 < L ; ix1++){
    for(int iy1 = 0; iy1 < L ; iy1++){
        int i = pos2index(0, ix1,  iy1);
        int ci = find(i);
        for(int ix2 = 0 ; ix2 < L; ++ix2){
        for(int iy2 = 0 ; iy2 < L; ++iy2){
            int j = (pos2index(L-1, ix2, iy2));
            int cj = find(j);
            if(ci == cj){
                return 1.0;
            }
        }
        }
    }
    }
    return 0.0;
}

int count_max_cluster(void) {
      std::vector<int> size(N, 0);
      for (int i = 0; i < N; i++) {
        int ci = find(i);
        size[ci]++;
      }
      int max = *std::max_element(size.begin(), size.end());
      return max;
}

int count_max_cluster_id(void) {
      std::vector<int> size(N, 0);
      for (int i = 0; i < N; i++) {
        int ci = find(i);
        size[ci]++;
      }
      auto max_el = std::max_element(size.begin(), size.end());
      int max_id = std::distance(size.begin(), max_el);
      return max_id;
}


int count_active_sites(void) {
    int sum = std::accumulate(sites.begin(), sites.end(), 0);
  return sum;
}


int count_trees(void) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        if(sites[i] == 1 and parent[i] == i){
            sum++ ;
        }
    }
    return sum;
}

void put_string_of_time(std::string &starting_date){
    std::time_t t = std::time(nullptr);
    std::stringstream ss ;
    ss<< std::put_time(std::localtime(&t), "%c %Z");
    starting_date+=ss.str();
}

int main(int argc, char** argv){
    // getting and printing time
    std::string starting_date;
    put_string_of_time(starting_date);
    std::cout << "Execution start at " << starting_date << std::endl;
  
    // making command line options
    std::string  fname_input;
    double  G0 ;
    std::string  fout_prefix;

    po::options_description opt("This program analyze percolation cluster");
    opt.add_options()
      ("help,h" ,                                          "show help")
      ("finp"   ,       po::value<std::string>(),          "file name of input")
      ("G0"     ,po::value<double>()->default_value(5.0),  "threshold of G, blow which sites to be connected")
      ("fout_prefix",   po::value<std::string>(),          "fout_prefix, a number of files will be made with this prefix  ");

     // analyze argc and argv and results are stored in vm
    try{
      po::variables_map vm;
      store(parse_command_line(argc, argv, opt), vm);
      notify(vm);
  
      if(vm.count("help")){
        std::cout << opt << std::endl; // show help
        exit(1);
      }
      else if(!vm.count("finp")){
        std::cerr << "finp is mandatory " << std::endl;
        std::cerr << "exit!!" << std::endl;
        exit(1);
      }
      else if(!vm.count("fout_prefix")){
        std::cerr << vm.count("fout_prefix") << std::endl;
        std::cerr << "fout_prefix is mandatory " << std::endl;
        std::cerr << "exit!!" << std::endl;
        exit(1);
      }
      else
      {
        fname_input = vm["finp"].as<std::string>();
        G0 = vm["G0"].as<double>();
        fout_prefix = vm["fout_prefix"].as<std::string>();
        std::cout << "**** Input parameters ****"  << std::endl;
        std::cout << "input file name: " << fname_input << std::endl;
        std::cout << "fout_prefix: " << fout_prefix << std::endl;
        std::cout << "G0: " << G0 << std::endl;
        std::cout << "*************************\n"  << std::endl;
      }
    }
    catch (boost::bad_any_cast &e) {
        std::cout <<"something wrong and buggy happened!!"  << std::endl;
        std::cout <<"exit!!"  << std::endl;
        exit(1);
    }
    catch (std::exception  &e) {
        std::cout << e.what() << std::endl;
        std::cout <<"exit!!"  << std::endl;
        exit(1);
    }


    int nixyz[3];

    //open fval files
    std::string buf;
    std::ifstream ifs_fmap(fname_input.c_str());
    if(!ifs_fmap){ std::cerr << "can not open " << fname_input << std::endl; exit(1);}
  
    // first line is number of grids in each direction
    std::getline(ifs_fmap, buf);
    sscanf(buf.c_str(), "%d %d %d", &nixyz[0], &nixyz[1], &nixyz[2]);

    if(nixyz[0] != nixyz[1] or nixyz[0] != nixyz[2]){
        std::cerr << "ERROR: currently, n_x = n_y = n_z should hold. exit!" << std::endl;
        return  EXIT_FAILURE;
    }
    init(nixyz[0]);

  
    //data follows
    while(std::getline(ifs_fmap, buf)){
        if(DEBUG) std::cout << buf <<std::endl;
  
        int  ix, iy, iz ;
        double tmp;
        sscanf(buf.c_str(), "%d %d %d %lf", &ix, &iy, &iz, &tmp);
    
        if( ix  < 0 or nixyz[0] <= ix ){
          std::cout << "error in ix " << ix << std::endl;
        }
        if( iy  < 0 or nixyz[1] <= iy ){
          std::cout << "error in iy " << iy << std::endl;
        }
        if( iz  < 0 or nixyz[2] <= iz ){
          std::cout << "error in iz " << iz << std::endl;
        }
    
        if(tmp > G0){
            sites[pos2index(ix,iy,iz)] = 0 ;
        }
        else{
            sites[pos2index(ix,iy,iz)] = 1 ;
        }
  
    }
    ifs_fmap.close();


    const std::string fname_log_name = fout_prefix+".log";

    std::ofstream ofs_log(fname_log_name.c_str());
    if(!ofs_log){
      std::cerr << "can not open: " << fname_log_name << std::endl;
      exit(1);
    }
      ofs_log << "Excecution starts at " << starting_date << "\n" <<std::endl;
  
    ofs_log <<  "executed command should be like: \n";
    for(int i = 0; argv[i] != NULL; i++){
      ofs_log << boost::format("%s ")% argv[i];
    }


    for (int ix = 0; ix < L; ++ix){
        for (int iy = 0; iy < L; ++iy){
            for (int iz = 0; iz < L; ++iz){
                int i = pos2index(ix,iy,iz);
                if(ix+1 < L){
                    connect(i, pos2index(ix+1, iy, iz));
                }
                if(iy+1 < L){
                    connect(i, pos2index(ix, iy+1, iz));
                }
                if(iz+1 < L){
                    connect(i, pos2index(ix, iy, iz+1));
                }
            }
        }
    } 

    std::cout << "G0 corss_x cross_y cross_z max_cluster #activesite #tree #sites maxid " <<std::endl;
    std::cout << G0 << " " ;
    std::cout << crossing_probability_x() << " " ;
    std::cout << crossing_probability_y() << " " ;
    std::cout << crossing_probability_z() << " " ;
    std::cout << count_max_cluster() << " " ;
    std::cout << count_active_sites() << " " ;
    std::cout << count_trees() << " " ;
    std::cout << N << " " ;
    std::cout << count_max_cluster_id() << " " ;
    std::cout << std::endl;

    const std::string fname_cluster_id = fout_prefix+"_cluster_id.dat";
    std::ofstream ofs_cluster(fname_cluster_id.c_str());
    if(!ofs_cluster){
      std::cerr << "can not open: " << fname_log_name << std::endl;
      exit(1);
    }
}
