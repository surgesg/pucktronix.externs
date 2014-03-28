/* 
-------------------------------------------------------------------------
ao_gsr - audio oracle object
implements an audio oracle
allows for generation of new sequences
with query informed by a second, parallel, data stream (gsr data)

copyright 2014 greg surges

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------
*/

#include <flext.h>
#if !defined(FLEXT_VERSION) || (FLEXT_VERSION < 400)
#error You need at least flext version 0.4.0
#endif

#define FLEXT_DEBUG 1

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

class ao_gsr:
	public flext_base
{
	FLEXT_HEADER(ao_gsr,flext_base)
 
public:
	ao_gsr();

protected:
    void m_init();
    void m_add_state(int, t_atom*);
    void m_add_sensor_data(float);
    void m_query(int, t_atom*);
    void m_print_sfx();
    void m_print_trn();
    void m_print_rsfx();
    void m_print_lrs();

    void m_next_state();
    void m_set_p(float);
    void m_set_lrs(int);
    void m_set_t(float);
    void m_set_n(int);
    void m_set_query(float);
    void m_toggle_regions(int);
    void m_set_region(int, int);

    void m_train_thresh(float, float, float);

    void m_save(t_symbol*);
    void m_load(t_symbol*);
    void m_get_ir();
    void load_oracle();

    int _find_better(int, vector<float>);
    float _get_distance(vector<float>, vector<float>);

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
    FLEXT_CALLBACK_V(m_add_state)
    FLEXT_CALLBACK_V(m_query)
    FLEXT_CALLBACK_I(m_toggle_regions)
    FLEXT_CALLBACK_I(m_set_lrs)
    FLEXT_CALLBACK_I(m_set_n)
    FLEXT_CALLBACK_II(m_set_region)
    FLEXT_CALLBACK_FFF(m_train_thresh)
    FLEXT_CALLBACK_F(m_set_p)
    FLEXT_CALLBACK_F(m_set_t)
    FLEXT_CALLBACK_F(m_set_query)
    FLEXT_CALLBACK_F(m_add_sensor_data)

    static void setup(t_classid);

    // navigation functions
    bool transition_from_k(int, int);
    vector<int> gather_options(int);
    vector<int> limit_to_region(vector<int>, int);

    // helper functions
    vector2D _encode();
    vector<float> _cumsum(vector<int>);
    vector<float> _calculate_ir();
    void _init_oracle();
    void _add_state(vector<float> new_event);

    // oracle data structures
    vector2D trn;
    vector1D sfx;
    vector2D rsfx;
    vector1D lrs;
    vector<vector<float> > data;
    vector<float> sensor_data;

    // oracle params
    int n_states;
    int max_lrs;
    float thresh;

    // navigation params
    float p; // jump prob
    int lrs_min; // constraint on navigation
    int n; // nav pointer
    int prev_n;
    int r_start, r_end;
    bool regions_active;
    float query_thresh;

};

FLEXT_NEW("ao_gsr",ao_gsr)

ao_gsr::ao_gsr(){ 
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
    FLEXT_ADDMETHOD_(0, "add_state", m_add_state);
    FLEXT_ADDMETHOD_F(0, "add_sensor_data", m_add_sensor_data);
    FLEXT_ADDMETHOD_(0, "query", m_query);
    FLEXT_ADDMETHOD_FFF(0, "train", m_train_thresh);

    FLEXT_ADDMETHOD_I(0, "set_lrs", m_set_lrs);
    FLEXT_ADDMETHOD_I(0, "set_n", m_set_n);
    FLEXT_ADDMETHOD_F(0, "set_p", m_set_p);
    FLEXT_ADDMETHOD_F(0, "set_t", m_set_t);
    FLEXT_ADDMETHOD_F(0, "set_query_thresh", m_set_query);
    FLEXT_ADDMETHOD_II(0, "set_region", m_set_region);

    m_init();
    thresh = 0.1;
    query_thresh = 0.1;

    post("ao_gsr - factor oracle object");
    post("2014 - greg surges");
    post("crel @ ucsd");
} 

void ao_gsr::m_init(){
    // initialize oracle and parameters
    _init_oracle();
    data.clear();
    data.push_back(vector<float>(0));

    sensor_data.clear();
    sensor_data.push_back(0);

    p = 0.5;
    lrs_min = 0;
    n = 0;
    prev_n = 1;
    r_start = 0;
    r_end = 1;
    regions_active = false;
}

void ao_gsr::_init_oracle(){
    // initialize oracle data structures (except data vec)
    sfx.clear();
    sfx.push_back(-1);

    vector1D temp;
    // temp.push_back(0);
    trn.clear();
    trn.push_back(temp);

    vector1D temp2;
    // temp2.push_back(0);
    rsfx.clear();
    rsfx.push_back(temp2);

    lrs.clear();
    lrs.push_back(0);
    
    n_states = 1; // including state 0?
    max_lrs = 0;
}

void ao_gsr::m_set_n(int new_n){
    // set oracle navigation parameter
    if((new_n >= 0) and (new_n < n_states)){
        n = new_n;
    }
}

void ao_gsr::m_set_query(float new_t){
    // set a min LRS constraint
    if(new_t >= 0){
        query_thresh = new_t;
    }
    post("query threshold set to: %f", new_t);
}

void ao_gsr::m_set_lrs(int new_lrs){
    // set a min LRS constraint
    if(new_lrs >= 0){
        lrs_min = new_lrs;
    }
    post("lrs minimum set to: %d", lrs_min);
}


void ao_gsr::m_set_t(float new_t){
    // set oracle distance threshold
    if(new_t > 0){
        thresh = new_t;
    } else {
        post("value of %f for t out of range", new_t);
    }
}

void ao_gsr::m_set_p(float new_p){
    // set oracle continuity parameter
    if((new_p > 0) and (new_p <= 1.0)){
        p = new_p;
    } else {
        post("value %f for p out of range", new_p);
    }
}

void ao_gsr::m_toggle_regions(int active){
    // toggle region constraints
    if(active == 1){
        post("regions active");
        regions_active = true;
    } else {
        post("regions inactive");
        regions_active = false;
    }
}

void ao_gsr::m_set_region(int start, int end){
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

float ao_gsr::_get_distance(vector<float> event1, vector<float> event2){
    // get distance between two events
    float distance = 0;
    for(int i = 0; i < event1.size(); i++){
        // post("event1[i] = %f, event2[i] = %f", event1[i], event2[i]);
        distance += ((event1[i] - event2[i]) * (event1[i] - event2[i]));
    }
    distance = sqrt(distance);
    return distance;
}

void ao_gsr::m_add_sensor_data(float new_datum){
    sensor_data.push_back(new_datum);
}


void ao_gsr::_add_state(vector<float> new_event){
    // add a state to oracle

    // add new default frame
    sfx.push_back(0);
    rsfx.push_back(vector<int>(0));
    trn.push_back(vector<int>(0));
    lrs.push_back(0);

    n_states += 1;
    int i = n_states - 1;
    trn[i-1].push_back(i);

    int k = sfx[i-1];
    int pi_1 = i-1;
    
    // add forward links
    vector<int> I;
    vector<float> dvec;
    while(k >= 0){
        dvec.clear();
        I.clear();
        for(int j = 0; j < trn[k].size(); j++){
            dvec.push_back(_get_distance(new_event, data[trn[k][j]]));
            // post("distance: %f", dvec[j]);
            if(dvec[j] < thresh){
                //I.push_back(trn[k][j]);
                I.push_back(j);
            }
        }
        if(I.size() == 0){
            trn[k].push_back(i);
            pi_1 = k;
            k = sfx[k];
        } else {
            break;
        }
    }
    //
    // if backtrack ended before last state 
    if (k >= 0){
        // get index of lowest value in I
        int argmin;
        vector<int> dvec_i;
        for(int j = 0; j < I.size(); j++){
            dvec_i.push_back(dvec[j]);
        }


        int min_val = 10000000;
        for(int j = 0; j < dvec_i.size(); j++){
            if(dvec_i[j] < min_val){
                argmin = j;
                min_val = dvec[j];
            }
        }
        
        sfx[i] = trn[k][I[argmin]];
        int S_i = sfx[i];
    } else {
        sfx[i] = 0;
    }

    k = _find_better(i, data[i-lrs[i]]);
    if(k != 0){
        lrs[i] += 1;
        sfx[i] = k;
    }
    rsfx[sfx[i]].push_back(i);
    
    if(lrs[i] > max_lrs){
        max_lrs = lrs[i];
    }
    
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
}

void ao_gsr::m_add_state(int argc, t_atom *argv){
    // wrapper for _add_state()
    // add new event/state to oracle

    vector<float> new_event;
    for(int j = 0; j < argc; j++){
        float new_dim;
        new_dim = GetAFloat(argv[j]);

        new_event.push_back(new_dim);
    }

    data.push_back(new_event);

    _add_state(new_event);
    ToOutInt(0, n_states);
}


int ao_gsr::_find_better(int i, vector<float> event){
    // from cheng-i's python code
    sort(rsfx[i].begin(), rsfx[i].end());
    for(int j = 0; j < rsfx[i].size(); j++){
        if((lrs[j] == lrs[i]) and _get_distance(data[j - lrs[i]], event) < thresh){
            return j;
        }
    }
}

void ao_gsr::m_query(int argc, t_atom *argv){
    // argument format:
    // size of feature vec, feature vec, sensor data
    //
    float rand_val = float(rand() % 100) / 100;
    // get dimensions of feature vector
    // this would be easier to just assume that the sensor data
    // is the last N atoms, and the feature is the atoms up to
    // that point...
    int feature_dim = GetAInt(argv[0]);
    // sensor input will be last
    float sensor_input = GetAFloat(argv[argc-1]);

    int new_n = n;
    if(rand_val > p){
        vector<float> new_event;
        for(int j = 1; j < feature_dim + 1; j++){
            float new_dim;
            new_dim = GetAFloat(argv[j]);
            new_event.push_back(new_dim);
        }

        int max_pos = 20;
        vector<int> connected_states;
        int k = n;
        // collect all states connected to current state
        while(k > 0){
            // add suffix
            connected_states.push_back(sfx[k]+1);
            // add rev suffixes
            for(int j = 0; j < rsfx[k].size(); j++){
                connected_states.push_back(rsfx[k][j] + 1);
                // this needs to sort of recusively add reverse suffixes as well
                for(int x = 0; x < rsfx[rsfx[k][j]].size(); x++){
                    connected_states.push_back(rsfx[rsfx[k][j]][x] + 1);
                }
            }
            // add trns
            for(int j = 0; j < trn[k].size(); j++){
                connected_states.push_back(trn[k][j] + 1);
            }
            // jump back along suffix
            k = sfx[k];
        }
        // remove duplicates
        sort(connected_states.begin(), connected_states.end());
        connected_states.erase(unique(connected_states.begin(),
                                      connected_states.end()),
                               connected_states.end());

        // collect states which are connected, below thresh, and in oracle bounds
        // in this case, we're checking for gsr similarity
        vector<int> intersection; 
        for(int i = 0; i < connected_states.size(); i++){
            // put sensor data into vectors
            vector<float> gsr_one, gsr_two;
            gsr_one.push_back(sensor_data[connected_states[i]]);
            gsr_two.push_back(sensor_input);

            if((_get_distance(gsr_one, gsr_two) < query_thresh) 
               && (connected_states[i] < n_states - 1)){
                intersection.push_back(connected_states[i]);
            }    
        }
        // choose closest matching state (maybe revise this?)
        // again, considering grs data as our target
        int nearest_index;
        float nearest_distance = 100000;
        float current_distance = 0;
        for(int i = 0; i < intersection.size(); i++){
            vector<float> gsr_one, gsr_two;
            gsr_one.push_back(sensor_data[intersection[i]]);
            gsr_two.push_back(sensor_input);
            current_distance = _get_distance(gsr_one, gsr_two); 
            if(current_distance < nearest_distance){
                nearest_distance = current_distance;
                nearest_index = i;
            }
        }
        int max_iterations = 10;
        int num_iterations = 0;
        if(intersection.size() > 0){
            while((new_n == n) and (num_iterations < max_iterations)){
                num_iterations += 1;
                new_n = intersection[nearest_index];
                if(n == 0){
                    new_n = intersection[rand() % intersection.size()];
                }
            }
            if(num_iterations == max_iterations){
                new_n = n+1;
            }
        } else {
            // random step or jump
            if((rand() % 100) > 50){
                new_n = n+1;
            } else {
                new_n = sfx[n];
            }
        }
    } else {
        new_n = n+1;
    }
    // if past bounds of oracle
    // choose a new random state to begin from
    if((new_n >= n_states) or (new_n == 0)){
        new_n = rand() % n_states;
    }
    n = new_n;
    ToOutInt(1, n);
}

void ao_gsr::m_next_state(){
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
    //ToOutInt(2, data[new_n]); // may want to add this in at some point
}

vector<int> ao_gsr::limit_to_region(vector<int> options, int n){
    // make new vector containing only those states which
    // fall into the selected region
    int max_iterations = 10;
    vector<int> new_options;
    
    for(int i = 0; i < options.size(); i++){
        if((options[i] >= r_start) and (options[i] <= r_end)){
            new_options.push_back(options[i]);
        }
    }
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

vector<int> ao_gsr::gather_options(int n){
    // gather possible jump locations and return them in a vector
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

void ao_gsr::m_train_thresh(float low_r, float high_r, float step){
    // find ideal threshold for audio oracle based on maximizing IR
    int num = (high_r - low_r) / step; // if n rounds down this could cause crash 
    float irs[num];
    double thresholds[num];
    int i = 0;
    vector<float> ir;
    vector<vector<float> > old_data = data;
    for(float t = low_r; t < high_r; t += step){
        m_init();

        thresholds[i] = t;
        thresh = t;

        for(int d = 1; d < old_data.size(); d++){
            _add_state(old_data[d]);
            data.push_back(old_data[d]);
        }
        ir = _calculate_ir();
        float ir_sum = 0;
        for(int s = 0; s < ir.size(); s++){
            ir_sum += ir[s];
        }
        irs[i] = ir_sum;
        i++;
    }
    // index of max ir val
    int max_ir_i = -1;
    int max_ir = -100000000;
    for(int i = 0; i < num; i++){
        if(irs[i] >= max_ir){
            max_ir = irs[i];
            max_ir_i = i;
        }
    }
    // m_init();
    // set thresh here
    post("ideal threshold set to %f", thresholds[max_ir_i]);
    thresh = thresholds[max_ir_i];

}

vector<float> ao_gsr::_calculate_ir(){
    // calculate IR
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
    vector<float> IR(n_states, 0);
    for(int i = 0; i < H0.size(); i++){
        H0[i] = log2(H0[i] + 0.0000001);
        H1[i] = log2(H1[i] + 0.0000001);
        H1[i] = H1[i] / BL[i];
        IR[i] = H0[i] - H1[i];
        // post("H0: %f, H1 %f, BL: %d, IR: %f", H0[i], H1[i], BL[i], IR[i]);
        if(IR[i] < 0){
            IR[i] = 0;
        }
    }
    return IR;
}

void ao_gsr::m_get_ir(){
    // calculate and output IR of oracle
    vector<float> IR = _calculate_ir();
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

vector<float> ao_gsr::_cumsum(vector<int> c0){
    // calculate cumulative sum of vector elements
    vector<float> res;
    int sum = 0;
    for(int i = 0; i < c0.size(); i++){
        sum += c0[i]; 
        res.push_back(sum);
    }
    return res;
}

vector2D ao_gsr::_encode(){
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

void ao_gsr::m_load(t_symbol *s){
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
    vector<vector<float> >new_data_vec;
    for(int i = 0; i < n_states; i++){
        string buf;
        vector1D new_data;
        getline(file, line);
        stringstream ss(line);
                                            
        vector<float> new_data_row;
        while(ss >> buf){
            float temp;
            stringstream(buf) >> temp;
            new_data_row.push_back(temp);
        }
        new_data_vec.push_back(new_data_row);
    }
    data = new_data_vec;
}

void ao_gsr::m_save(t_symbol *s){
    // save an oracle to a text file
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
        for(int j = 0; j < data[i].size(); j++){
            file << data[i][j] << endl;
        }
    }
    file.close();
    post("done");
}


void ao_gsr::m_print_sfx(){
    // print suffixes
    post("dumping suffixes:");
    for(int i = 0; i < sfx.size(); i++){
        post("state %d: %d", i, sfx[i]);
    }

}

void ao_gsr::m_print_trn(){
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

void ao_gsr::m_print_rsfx(){
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

void ao_gsr::m_print_lrs(){
    // print LRS vector
    post("dumping lrs:");
    for(int i = 0; i < lrs.size(); i++){
        post("state %d: %d", i, lrs[i]);
    }

}





