/* 
-------------------------------------------------------------------------
fo - factor oracle object
implements a symbolic factor oracle
allows for generation of new sequences

copyright 2014 greg surges

-------------------------------------------------------------------------
*/

#include <flext.h>
#if !defined(FLEXT_VERSION) || (FLEXT_VERSION < 400)
#error You need at least flext version 0.4.0
#endif

//#define FLEXT_DEBUG 1

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;
typedef vector<int> vector1D;
typedef vector<vector1D> vector2D;

class fo:
	public flext_base
{
	FLEXT_HEADER(fo,flext_base)
 
public:
	fo();

protected:
    void m_init();
    void m_add_state(int);
    void m_print_sfx();
    void m_print_trn();
    void m_print_rsfx();
    void m_print_lrs();

    void m_next_state();
    void m_set_p(float);
    void m_set_lrs(int);
    void m_set_n(int);
    void m_toggle_regions(int);
    void m_set_region(int, int);

    void m_save(t_symbol*);
    void m_load(t_symbol*);
    vector<float> m_get_ir();
    void load_oracle();

    int _find_better(int, int);

private:
    FLEXT_CALLBACK(m_init)
    FLEXT_CALLBACK(m_print_sfx)
    FLEXT_CALLBACK(m_print_trn)
    FLEXT_CALLBACK(m_print_rsfx)
    FLEXT_CALLBACK(m_print_lrs)
    FLEXT_CALLBACK(m_next_state)
    FLEXT_CALLBACK(m_get_ir)
    FLEXT_CALLBACK_S(m_save)
    FLEXT_CALLBACK_S(m_load)
    FLEXT_CALLBACK_I(m_add_state)
    FLEXT_CALLBACK_I(m_set_lrs)
    FLEXT_CALLBACK_I(m_set_n)
    FLEXT_CALLBACK_I(m_toggle_regions)
    FLEXT_CALLBACK_II(m_set_region)
    FLEXT_CALLBACK_F(m_set_p)

    static void setup(t_classid);

    // navigation functions
    bool transition_from_k(int, int);
    vector<int> gather_options(int);
    vector<int> limit_to_region(vector<int>, int);

    // helper functions
    vector2D _encode();
    vector<float> _cumsum(vector<int>);

    // oracle data structures
    vector2D trn;
    vector1D sfx;
    vector2D rsfx;
    vector1D lrs;
    vector1D data;

    // oracle params
    int n_states;
    int max_lrs;

    // navigation params
    float p; // jump prob
    int lrs_min; // constraint on navigation
    int n; // nav pointer
    int r_start, r_end;
    bool regions_active;
};

FLEXT_NEW("fo",fo)

fo::fo(){ 
    // init oracle
	AddInAnything();  // add one inlet for any message
	
	// define outlets:
	AddOutInt();  // add one int outlet (has index 0)
	AddOutInt();  // add one int outlet (has index 1)
	AddOutInt();  // add one int outlet (has index 2)
	AddOutAnything();  // add one int outlet (has index 3)
	
	// register methods
    FLEXT_ADDMETHOD_(0, "init", m_init);

    FLEXT_ADDMETHOD_(0, "print_sfx", m_print_sfx);
    FLEXT_ADDMETHOD_(0, "print_trn", m_print_trn);
    FLEXT_ADDMETHOD_(0, "print_lrs", m_print_lrs);
    FLEXT_ADDMETHOD_(0, "print_rsfx", m_print_rsfx);

    FLEXT_ADDMETHOD_(0, "next_state", m_next_state);
    FLEXT_ADDMETHOD_(0, "get_ir", m_get_ir);

    FLEXT_ADDMETHOD_(0, "save", m_save);
    FLEXT_ADDMETHOD_(0, "load", m_load);

    FLEXT_ADDMETHOD_I(0, "toggle_regions", m_toggle_regions);
    FLEXT_ADDMETHOD_I(0, "add_state", m_add_state);

    FLEXT_ADDMETHOD_I(0, "set_lrs", m_set_lrs);
    FLEXT_ADDMETHOD_F(0, "set_p", m_set_p);
    FLEXT_ADDMETHOD_F(0, "set_n", m_set_n);
    FLEXT_ADDMETHOD_II(0, "set_region", m_set_region);

    post("fo - factor oracle object");
    post("2014 - greg surges");
    post("crel @ ucsd");
    
    m_init();
} 

void fo::m_init(){
    // initialize oracle and parameters
    sfx.clear();
    sfx.push_back(-1);

    vector1D temp;
    trn.clear();
    trn.push_back(temp);

    vector1D temp2;
    rsfx.clear();
    rsfx.push_back(temp2);

    lrs.clear();
    lrs.push_back(0);
    data.clear();
    data.push_back(0);

    n_states = 1;
    max_lrs = 0;

    p = 0.5;
    lrs_min = 0;
    n = 0;
    r_start = 0;
    r_end = 1;
    regions_active = false;
}

void fo::m_set_n(int new_n){
    // set oracle navigation parameter
    if((new_n >= 0) and (new_n < n_states)){
        n = new_n;
    }
}

void fo::m_set_p(float new_p){
    // set oracle continuity parameter
    if((new_p > 0) and (new_p <= 1.0)){
        p = new_p;
    } else {
        post("value %f for p out of range", new_p);
    }
}

void fo::m_set_lrs(int new_lrs){
    // set a min LRS constraint
    if(new_lrs >= 0){
        lrs_min = new_lrs;
    }
    post("lrs minimum set to: %d", lrs_min);
}

void fo::m_toggle_regions(int active){
    // toggle region constraints
    if(active == 1){
        post("regions active");
        regions_active = true;
    } else {
        post("regions inactive");
        regions_active = false;
    }
}

void fo::m_set_region(int start, int end){
    // set region constraints
    if(start < 0){
        start = 0;
    }
    if(end >= n_states){
        end = n_states - 1;
    }
    if(start >= end){
        end = start + 1;
    }
    r_start = start;
    r_end = end;
    post("region set to (start: %d, end: %d)", start, end);
}

void fo::m_add_state(int new_symbol){
    // add a state to oracle

    // add new default frame
    sfx.push_back(0);
    data.push_back(new_symbol);
    rsfx.push_back(vector<int>(0));
    trn.push_back(vector<int>(0));
    lrs.push_back(0);

    n_states += 1;
    int i = n_states - 1;
    trn[i-1].push_back(i);

    int k = sfx[i-1];
    int pi_1 = i-1;
    
    // add forward links
    while((k >= 0) and !transition_from_k(k, new_symbol)){
        trn[k].push_back(i);
        pi_1 = k;
        k = sfx[k];
    }

    // if backtrack ended before last state 
    if (k >= 0){
        int S_i = 0;
        for(int t = 0; t < trn[k].size(); t++){
            if(data[trn[k][t]] == new_symbol){
               S_i = trn[k][t]; 
            }
        }
        sfx[i] = S_i;
    } else {
        sfx[i] = 0;
    }

    k = _find_better(i, data[i-lrs[i]]);
    if(k != 0){
        sfx[i] = k;
    }
    rsfx[sfx[i]].push_back(i);

    // LRS
    int ss = sfx[i];
    if((ss == 0) or (ss == 1)){
        lrs[i] = 0;
    } else {
        int pi_2 = ss - 1;
        if(pi_2 == sfx[pi_1]){
            lrs[i] = lrs[pi_1] + 1;
        } else {
            while(sfx[pi_2] != sfx[pi_1]){
                pi_2 = sfx[pi_2];
            }
            if(lrs[pi_1] < lrs[pi_2]){
                lrs[i] = lrs[pi_1] + 1;
            } else {
                lrs[i] = lrs[pi_2] + 1;
            }
        }
    }
    ToOutInt(0, n_states);
}

bool fo::transition_from_k(int k, int new_symbol){
    // check if k has a transition labeled by new_symbol
    bool trnpresent = false;
    vector<int> symbols = trn[k];
    // check for transition to new data from k
    for(int j = 0; j < symbols.size(); j++){
        if(new_symbol == data[symbols[j]]){
            trnpresent = true;
        }
    }
    return trnpresent;
}

int fo::_find_better(int i, int symbol){
    // from cheng-i's python code
    sort(rsfx[i].begin(), rsfx[i].end());
    for(int j = 0; j < rsfx[i].size(); j++){
        if((lrs[j] == lrs[i]) and (data[j - lrs[i]] == symbol)){
            return j;
        }
    }
}

void fo::m_next_state(){
    // make one step of navigation along oracle
    int new_n;
    float rand_val = float(rand() % 1000) / 1000;
    if(n_states <= 1){
        post("invalid oracle");
        return;
    }
    if((rand_val < p) and (n < n_states - 1)){
        // transition
        new_n = n + 1;
        if(regions_active){
            if((new_n < r_start) or (new_n > r_end)){
                new_n = (rand() % (r_end - r_start)) + r_start;
            }
        }
    } else {
        // jump
        vector<int> options = gather_options(n);
        if(regions_active){
            // remove states outside region
            options = limit_to_region(options, n);
        }
        if(options.size() > 0){
            // choose new state, making sure it's not
            // an immediate repetition
            int rand_choice = rand() % options.size();
            new_n = options[rand_choice];
        } else {
            // if no connected states exist, do a transition to the
            // next state. however, if that state is out of region,
            // then reset to beginning of region
            if(new_n == r_end){
                new_n = r_start;
            } else {
                new_n = n + 1;
            }
        }
    }
    // REMOVE IMMEDIATE REPETITIONS
    if(new_n == n){
        new_n = n+1;
    }
    // BOUNDARY CONDITIONS
    if(new_n >= n_states){
        new_n = rand() % n_states;
    }
    if((new_n == 0) or (new_n == 1)){
        // choose a transition from 0
        if(regions_active){
            vector<int> r_options = trn[0];
            r_options = limit_to_region(r_options, n);
            new_n = r_options[rand() % r_options.size()];
        } else {
            new_n = trn[0][rand() % trn[0].size()];
        }
    }
    n = new_n;
    ToOutInt(1, new_n);
    ToOutInt(2, data[new_n]);
}

vector<int> fo::limit_to_region(vector<int> options, int n){
    // make new vector containing only those states which
    // fall into the selected region
    int max_iterations = 10;
    vector<int> new_options;

    for(int i = 0; i < options.size(); i++){
        if((options[i] >= r_start) and (options[i] <= r_end)){
            new_options.push_back(options[i]);
        }
    }
    //
    // if none remain, choose one at random
    if(new_options.size() == 0){
        int new_val = n;
        int n_iterations = 0;
        while((new_val == n) and (n_iterations < max_iterations)){ 
            new_val = (rand() % (r_end - r_start)) + r_start;
            n_iterations += 1;
        }
        if(n_iterations >= max_iterations){
            new_val = r_start;
        }
        new_options.push_back(new_val);
    }
    return new_options;
}

vector<int> fo::gather_options(int n){
    // gather possible jump locations from n and return them in a vector
    // add suffix
    vector<int> options;
    options.push_back(sfx[n] + 1); // suffix

    // add reverse suffix
    int rsfx_size = rsfx[n].size(); // hold this steady in case modified elsewhere
    for(int j = 0; j < rsfx_size; j++){
        if(rsfx[n][j] != n+1){ // if not redundant
            options.push_back(rsfx[n][j] + 1); 
        }
    }

    // suffix of suffix
    options.push_back(sfx[sfx[n]] + 1);

    // rsfx of rsfx
    for(int j = 0; j < rsfx_size; j++){ // all rev sfx of current
        for(int r = 0; r < rsfx[rsfx[n][j]].size(); r++){ // rev sfx of rev sfx
            int cur_val = rsfx[rsfx[n][j]][r];
            if(cur_val != n + 1){
                options.push_back(cur_val);
            }
        }
    }
    // remove duplicates
    sort(options.begin(), options.end());
    options.erase(unique(options.begin(), options.end()), options.end());
    // filter for n+1 and LRS
    for(int j = 0; j < options.size(); j++){
        if((options[j] == n+1) or (lrs[options[j]] < lrs_min)){
            options.erase(options.begin() + j);
        }
    }
    return options;
}

vector<float> fo::m_get_ir(){
    // calculate IR of oracle and output
    vector2D code = _encode();
    vector<int> cw0(n_states, 0);
    vector<int> cw1(n_states, 0);
    vector<int> BL(n_states, 0);

    int j = 0;
    for(int i = 0; i < code.size(); i++){
        if(code[i][0] == 0){
            cw0[j] = 1;
            cw1[j] = 1;
            BL[j] = 1;
            j += 1;
        } else {
            int L = code[i][0];
            for(int k = j; k < j+L; k++){
                cw0[k] = 0;
                if(k == j){
                    cw1[k] = 1;
                } else {
                    cw1[k] = 0;
                }
                BL[k] = L;
            }
            j += L;
        }
    }

    vector<float> H0 = _cumsum(cw0);
    vector<float> H1 = _cumsum(cw1);
    float IR[n_states];
    for(int i = 0; i < H0.size(); i++){
        H0[i] = log2(H0[i] + 0.0000001);
        H1[i] = log2(H1[i] + 0.0000001);
        H1[i] = H1[i] / BL[i];
        IR[i] = H0[i] - H1[i];
        if(IR[i] < 0){
            IR[i] = 0;
        }
    }

    AtomList result(n_states+1); // one for "list" argument
    SetSymbol(result[0], sym_list);
    t_atom new_atom;
    for(int i = 1; i < n_states; i++){
        SetFloat(new_atom, IR[i]);
        CopyAtom(&result[i+1], &new_atom); // offset to compensate for "list" argument
    }
    // output result list as an anything
    ToOutAnything(3, GetSymbol(result[0]), n_states, result.Atoms()+1);
}

vector<float> fo::_cumsum(vector<int> c0){
    vector<float> res;
    int sum = 0;
    for(int i = 0; i < c0.size(); i++){
        sum += c0[i]; 
        res.push_back(sum);
    }
    return res;
}

vector2D fo::_encode(){
    // encode oracle using compror algorithm and return compressed oracle
    vector2D compror;
    vector2D code;
    
    int j = 0;
    int i = j;
    int cnt = 1;
    while(j < lrs.size() - 1){
        while((i < lrs.size() - 1) and (lrs[i+1] >= i - j + 1)){
            i += 1;
        }
        if(i == j){
            i += 1;
            vector<int> new_code;
            new_code.push_back(0);
            new_code.push_back(i);
            code.push_back(new_code);
        } else {
            vector<int> new_code;
            new_code.push_back(i-j);
            new_code.push_back(sfx[i] - i + j + 1);
            code.push_back(new_code);

            vector<int> new_compror;
            new_compror.push_back(i);
            new_compror.push_back(i-j);
            compror.push_back(new_compror);
        }
        cnt += 1;
        j = i;
    }
    return code;
}

void fo::m_load(t_symbol *s){
    // load an oracle from .txt file
    string fname = GetString(s);

    post("loading %s", fname.c_str());
    ifstream file;
    file.open(fname.c_str());
    string line;
    getline(file, line);

    stringstream(line) >> n_states;
    post("loading %d states...", n_states);

    // first load sfxs
    vector1D new_sfx_vec;
    int new_sfx = 0;
    for(int i = 0; i < n_states; i++){
        getline(file, line);
        stringstream(line) >> new_sfx;
        new_sfx_vec.push_back(new_sfx);
    }
    sfx = new_sfx_vec;
    
    // then transitions
    vector2D new_trn_vec;
    for(int i = 0; i < n_states; i++){
        string buf;
        vector1D new_trn;
        getline(file, line);
        stringstream ss(line);

        vector1D new_trn_row;
        while(ss >> buf){
            int temp;
            stringstream(buf) >> temp;
            new_trn_row.push_back(temp);
        }
        new_trn_vec.push_back(new_trn_row);
    }
    trn = new_trn_vec;
    //
    // then rsfx 
    vector2D new_rsfx_vec;
    for(int i = 0; i < n_states; i++){
        string buf;
        vector1D new_rsfx;
        getline(file, line);
        stringstream ss(line);
                                            
        vector1D new_rsfx_row;
        while(ss >> buf){
            int temp;
            stringstream(buf) >> temp;
            new_rsfx_row.push_back(temp);
        }
        new_rsfx_vec.push_back(new_rsfx_row);
    }
    rsfx = new_rsfx_vec;
   
    // then lrs
    vector1D new_lrs_vec;
    int new_lrs = 0;
    for(int i = 0; i < n_states; i++){
        getline(file, line);
        stringstream(line) >> new_lrs;
        new_lrs_vec.push_back(new_lrs);
    }
    lrs = new_lrs_vec;

    // then data 
    vector1D new_data_vec;
    int new_data = 0;
    for(int i = 0; i < n_states; i++){
        getline(file, line);
        stringstream(line) >> new_data;
        new_data_vec.push_back(new_data);
    }
    data = new_data_vec;
}

void fo::m_save(t_symbol *s){
    // save an oracle to a txt file
    string fname = GetString(s);

    post("saving...");
    ofstream file;
    file.open(fname.c_str());
    // first line contains length of oracle
    file << n_states << endl;
    // then suffixes
    for(int i = 0; i < sfx.size(); i++){
        file << sfx[i] << endl;
    }
    // then transitions - each state on one line, divided by spaces
    for(int i = 0; i < trn.size(); i++){
        for(int j = 0; j < trn[i].size(); j++){
            file << trn[i][j] << " ";
        }
        file << endl;
    }
    // then rsfxs - each state on one line, divided by spaces
    for(int i = 0; i < rsfx.size(); i++){
        for(int j = 0; j < rsfx[i].size(); j++){
            file << rsfx[i][j] << " ";
        }
        file << endl;
    }
    // then lrs - each state on one line
    for(int i = 0; i < lrs.size(); i++){
        file << lrs[i] << endl;
    }
    // then data - each state on one line
    for(int i = 0; i < data.size(); i++){
        file << data[i] << endl;
    }
    file.close();
    post("done");
}


void fo::m_print_sfx(){
    // print suffixes
    post("dumping suffixes:");
    for(int i = 0; i < sfx.size(); i++){
        post("state %d: %d", i, sfx[i]);
    }

}

void fo::m_print_trn(){
    // print transitions
    post("dumping transitions:");
    for(int i = 0; i < trn.size(); i++){
        string composite;
        string new_str;
        for(int j = 0; j < trn[i].size(); j++){
            stringstream ss;
            ss << trn[i][j];
            new_str = ss.str();
            new_str += ", ";
            composite += new_str;
        }
        post("state %d: %s", i, composite.c_str());
    }

}

void fo::m_print_rsfx(){
    // print reverse suffixes
    post("dumping reverse suffixes:");
    for(int i = 0; i < rsfx.size(); i++){
        string composite;
        string new_str;
        for(int j = 0; j < rsfx[i].size(); j++){
            stringstream ss;
            ss << rsfx[i][j];
            new_str = ss.str();
            new_str += ", ";
            composite += new_str;
        }
        post("state %d: %s", i, composite.c_str());
    }
}

void fo::m_print_lrs(){
    // print lrs vector
    post("dumping lrs:");
    for(int i = 0; i < lrs.size(); i++){
        post("state %d: %d", i, lrs[i]);
    }

}






