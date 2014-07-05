// LMS-spec~ - a flext tutorial external written by Frank Barknecht
// 
// This is a commented port of the pan~ example from the PD-Externals-Howto to
// illustrate the usage of flext. You can get the original code at
// http://iem.kug.ac.at/pd/externals-HOWTO/

#include <flext.h>
#include <stdlib.h>
 
#if !defined(FLEXT_VERSION) || (FLEXT_VERSION < 401)
#error You need at least flext version 0.4.1
#endif


// A flext dsp external ("tilde object") inherits from the class flext_dsp 
class LMSspec: 
	public flext_dsp
{
	// Each external that is written in C++ needs to use #defines 
	// from flbase.h
	// 
	// The define
	// 
	// FLEXT_HEADER(NEW_CLASS, PARENT_CLASS)
	// 
	// should be somewhere in your dsp file.
	// A good place is here:
	
	FLEXT_HEADER(LMSspec, flext_dsp)

	public:
		LMSspec(int argc, t_atom *argv)
		{
            // arguments [LMSspec~ u L]
            // or just [LMSspec~ u] 

			AddInSignal("filter audio in");       // left audio in
			AddInSignal("source audio in");      // right audio in
			AddInAnything("messages");    // 1 float in
			AddOutSignal("audio out");          // 1 audio out 
			AddOutSignal("audio out");          // 1 audio out 
			AddOutAnything("data out");          // 1 audio out 
			
			// Now we need to bind the handler function to our
			// panning inlet, which is inlet 2 (counting all inlets
			// from 0).  We want the function "setPan" to get
			// called on incoming float messages:

			FLEXT_ADDMETHOD_(2, "printw", m_printw);
			FLEXT_ADDMETHOD_(2, "init", m_init);
            L = 16;
            u = 0.01;// a / (r0*L);

            switch(argc){
                case 2:
                    L = GetInt(argv[1]);
                    u = GetFloat(argv[0]);
                case 1:
                    u = GetFloat(argv[0]);
                case 0:
                    ;
            }

            w = new double[L];
            taps = new double[L];
            y = new double[L];
            for(int i = 0; i < L; i++){
                w[i] = taps[i] = y[i] = 0.f;
            }

            src_taps = new double[L];
            src_w = new double[L+1];
            for(int i = 0; i < L; i++){
                src_taps[i] = src_w[i] = 0.f;
            }
            src_w[L] = 0.f;
		
			// We're done constructing:
			post("-- LMSspec~ -- \n 2014 greg surges \n surgesg@gmail.com");
            post("L: %d, u: %f", L, u);
			
		} // end of constructor
		
	
	protected:
		// here we declare the virtual DSP function
		virtual void m_signal(int n, float *const *in, float *const *out);
        void m_printw();
        void m_init();
        ~LMSspec(); // destructor
	private:	
        float u;        // learning rate
        int L;          // num taps
        double *w;       // weight vector
        double *taps;    // delay line
        double *y;       // filter state

        double *src_taps;
        double *src_w;
        
		// Before we can use "setPan" as a handler, we must register this
		// function as a callback to PD or Max. This is done using the
		// FLEXT_CALLBACK* macros. There are several of them.
		//
		// FLEXT_CALLBACK_F is a shortcut, that registers a function
		// expecting one float arg (thus ending in "_F"). There are
		// other shortcuts that register other types of functions. Look
		// into flext.h. No semicolon at the end of line!!!
		FLEXT_CALLBACK(m_printw)
		FLEXT_CALLBACK(m_init)

}; // end of class declaration for LMS-spec


// Before we can run our LMS-spec-class in PD, the object has to be registered as a
// PD object. Otherwise it would be a simple C++-class, and what good would
// that be for?  Registering is made easy with the FLEXT_NEW_* macros defined
// in flext.h. For tilde objects without arguments call:

FLEXT_NEW_DSP_V("LMSspec~ LMSspec", LMSspec)
// T.Grill: there are two names for the object: LMS-spec~ as main name and pan~ as its alias

// Now we define our DSP function. It gets this arguments:
// 
// int n: length of signal vector. Loop over this for your signal processing.
// float *const *in, float *const *out: 
//          These are arrays of the signals in the objects signal inlets rsp.
//          oulets. We come to that later inside the function.

LMSspec::~LMSspec(){
    delete[] w;
    delete[] taps;
    delete[] y;

}

void LMSspec::m_init(){
    for(int i = 0; i < L; i++){
        w[i] = taps[i] = y[i] = 0.f;
    }
}

void LMSspec::m_printw(){
    // this outputs coefficients for 
    // instantaneous spectral estimate
    AtomList out(L+2);
    SetSymbol(out[0], sym_list);
    SetFloat(out[1], 1.0);
    t_atom new_atom;
    for(int i = 0; i < L; i++){
        //post("%d: %f", i, w[i]);
        SetFloat(new_atom, w[i]*-1);
        CopyAtom(&out[i+2], &new_atom);
    }
    ToOutAnything(2, GetSymbol(out[0]), L+1, out.Atoms()+1);
}

void LMSspec::m_signal(int n, float *const *in, float *const *out)
{
	
	//const float *ins1    =  in[0];
	const float *ins1    =  InSig(0);
	//const float *ins2    =  in[1];
	const float *ins2    =  InSig(1);
	// As said above "in" holds a list of the signal vectors in all inlets.
	// After these two lines, ins1 holds the signal vector ofthe first
	// inlet, index 0, and ins2 holds the signal vector of the second
	// inlet, with index 1.
	
	float *outs1          = OutSig(0);
	float *outs2          = OutSig(1);
	// Now outs holds the signal vector at the one signal outlet we have.
	
	// We are now ready for the main signal loop
	while (n--)
	{
        double dh = 0;   // estimate
        double e = 0;    // error
        double x = (*ins1++);

        // filter - dot product of weights and taps
        for(int i = 0; i < L; i++){
            dh += (taps[i] * w[i]); 
        }
        // get error
        e = x - dh;

        // update weights 
        // seems like it might be working, but as the input changes,
        // the spectral envelope doesn't change correctly
        // i.e. if a single sinusoid is input, if it changes frequency
        // the zero at the original freq will remain
        // this causes the error signal to remain zero/small and therefore
        // the false peak never adapts away
        float update[L];
        float rand_term;
        for(int i = 0; i < L; i++){
            //rand_term = (rand() % 10000) - 5000;
            //rand_term = float(rand_term) / 100000000; // 0 - 0.001
            // this code can be simplified
            update[i] = 2.0*u*e*taps[i];
            w[i] = w[i] + 2.0*u*e*taps[i]; // + rand_term;
        }
	    // update delay line
        for(int i = L-1; i > 0; i--){
            taps[i] = taps[i-1];
        }
        taps[0] = x; 
        // error sig
		*outs1++  = e;

        //**** begin source/filter method here *****************
        double s = (*ins2++);
        double src_filtered = 0;

        // process zeros to poles / FIR to IIR
        src_w[0] = 1;
        for(int i = 0; i < L; i++){
            src_w[i+1] = -w[i];
        }
        
        src_filtered = -1 * s; // B
        for(int i = 1; i < L+1; i++){ // As
            src_filtered -= src_w[i] * src_taps[i]; 
        }

        for(int i = L; i > 0; i--){
            src_taps[i] = src_taps[i-1];
        }
        src_taps[0] = src_filtered; 

        *outs2++ = src_filtered;
        //*outs2++ = s;
	}
}  // end m_signal
