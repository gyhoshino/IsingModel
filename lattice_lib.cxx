// possible improvements:
// calculate the common values of deltaE (depending on number of bonds) once instead of repeatedly throughout the run
// float deltaE = 2*(J*site[*r].n_bonds+B*(2*site[*r].spin-1));
// if (deltaE < 0.0 || dis(gen) < exp(-1.0*deltaE/T))
//
// compile with optimization -O1 or -O2 or -O3
//
// use the Wolff algorithm with a recursive implementation
// (much of the implementation taken from:
// https://www.uio.no/studier/emner/matnat/fys/nedlagte-emner/FYS4410/v06/undervisningsmateriale/Programs/Week\%2013/wolff.cpp

#include <random>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
 
using namespace std;

struct lattice_site{
    /* ~lattice_site(){ */ 
        /* cout << "del lattice" << endl; */
        /* delete left; delete right; delete up; delete down; }; */
    /* ~lattice_site(){ delete left; delete right; delete up; delete down; }; */
	bool  spin;       // whether it is spin up (true) or spin down (false)
	short n_bonds;    // number of aligned spins with neighbor (+1 aligned, -1 anti-aligned)
	bool  bond_right; // true if spin aligned with site to right
	bool  bond_down;  // true is spin aligned with site down
    /* int */  
    lattice_site  *left;
    lattice_site  *right;
    lattice_site  *up;
	lattice_site  *down;
};

class IsingMatrix{
    // Make it of lattice sites (ls), each of which is linked to four other 
    // sites: up/down/left/rigth: U/D/L/R.

    // Each site "owns" the boolean of opposing flip with the site next to
    // it on the D/R sides

	// Each LS has a pointer to each four adjacent sites


	// Each site also has pointers to the each of the sites 
    public:
    IsingMatrix(int set_N = 10, int set_Npick = 10, int offset = 0); // same as allocate
    ~IsingMatrix();
    
    float get_E() { return -(J*Nbond+B*Nspin)/NN; };
    float get_M() { return (float)Nspin/NN; };
    int   get_N() { return N; };
    int   set_Npick(int in_val) { Npick = in_val; return 0; };
    int   set_J (float in_val) { J = in_val; return 0;} ;
    int   print_E() { cout << " Value of E: " << get_E() << endl; return 0; };
    int   print_Nspin() { cout << " #spin: " << Nspin << endl; return 0; };
    int   get_Nspin() { return Nspin; };
    int   get_Nbond() { return Nbond; };
    int   print_Nbond() { cout << " #bonds: " << Nbond << endl; return 0; };
    int   get_spin(int i, int j) { return (site[i*N + j].spin ? 1 : -1); };
    float get_site_flip_prob(int i, int j, float T, float B);
    inline int   get_spin(int i) { return (site[i].spin ? 1 : -1); };
    int   set_spin(int i, int j, int k) ;
    void  flip_site(lattice_site& _site);
          /* { site[i*N + j].spin = ((k > 0) ? true : false); return 0; }; */
    int   print_spins();     // defined
    int   print_bonds();     // defined
    int   initialize_spins(); // defined

    float T_step_prior; // initialize to -1, thereafter use T_step_prior adn B_step_prior to determine if
    float B_step_prior; // new flip_prob values need to be calculated
    float flip_prob[5][2]; // for any given T and B, there are only five distinct flip prob.
                        // swapping 4 bonds(x2:up/down), 2 bonds(x2:up/down), 0 bonds
    void  update_flip_prob(float T, float B);
    int   step(float T, float B);

    // Note: the results of using wolff are not expected. 
    // Do not use until it is debugged.
    int   wolff(float T); // evolve the lattice with a Wolff step
    float w_clust_prob;   // clustering probability for wolff cluster == add if dis(gen) < 1 - exp(-2/T) // same as dis(gen) > exp(-2/T)
    bool  w_cluster_spin;
    void  w_grow_cluster(lattice_site&);
    void  w_try_add_site(lattice_site*);

    int   nsteps(int nsteps, float T, float B);
    int   calc_auto_correlation();
    float *auto_correlation;
    int n_wolff;

    int   n_spincorr;
    int   calc_spincorr();
    float *spincorr;

    /* bool *spin; */
    /* short int *pairs; */

    private:
	lattice_site* site;

    int N;
    int NN;
    int Nbond;
    int Nspin;
    float J = 1.0;
    float E;
    float M;
    float B;
    /* std::random_device rd; */
    std::mt19937 gen;
    std::uniform_real_distribution<float> dis;

    vector<int> indices; // used to pick the order of lattice sights
    int Npick;
};

int IsingMatrix::wolff(float T){
    /* n_wolff = 0; */
    B = 0;
    /* update_flip_prob(T,B); */

    /* w_clust_prob = exp(-2/T); */
    int i_start = (int)(dis(gen)*NN);
    w_clust_prob = 1 - exp(-2/T);
    /* w_clust_prob = 1.; */
    /* cout << "n_wolff start " << i_start/N << " " << i_start %N<< endl; */
    // get a random starting site
    bool w_cluster_spin{site->spin};
    w_grow_cluster(site[i_start]);
    return 0;
};
void IsingMatrix::w_grow_cluster(lattice_site& _site){
    ++n_wolff;
    flip_site(_site);
    /* cout << "n_wolff " << n_wolff << endl; */
    /* print_spins(); */
    w_try_add_site(_site.left);
    w_try_add_site(_site.right);
    w_try_add_site(_site.up);
    w_try_add_site(_site.down);
};
void IsingMatrix::w_try_add_site(lattice_site* _site){
    if ( (_site->spin == w_cluster_spin)
      && (dis(gen) < w_clust_prob)
    ) w_grow_cluster(*_site);

};

float IsingMatrix::get_site_flip_prob(int i, int j, float T, float B) {
    int n = i*N+j;
    update_flip_prob(T,B);
    return flip_prob[2+site[n].n_bonds/2][site[n].spin];
};

int IsingMatrix::set_spin(int ii, int jj, int kk){
    int i = ii*N + jj;
    bool set_spin_state = kk > 0;
    if (site[i].spin != set_spin_state) flip_site(site[i]);
    return 1;
};

void IsingMatrix::flip_site(lattice_site& _site){
    Nbond -= 2*_site.n_bonds;
    Nspin += ( _site.spin ? -2 : +2 );
  
    _site.spin    ^= true;
    _site.n_bonds *= -1;
    
    _site.right->n_bonds += ( _site.bond_right ? -2 : 2);
    _site.bond_right     ^= true;
  
    _site.down->n_bonds += ( _site.bond_down ? -2 : 2);
    _site.bond_down     ^= true;
  
    _site.up->n_bonds   += ( _site.up->bond_down ? -2 : 2);
    _site.up->bond_down ^= true;
  
    _site.left->n_bonds    += ( _site.left->bond_right ? -2: 2);
    _site.left->bond_right ^= true;
}

int IsingMatrix::print_spins() {
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            cout << " " << site[i*N+j].spin;
        }
        cout << endl;
    };
    return 0;
}
int IsingMatrix::print_bonds() {
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            printf("%3i ", site[i*N+j].n_bonds);
            /* cout << " " << site[i*N+j].n_bonds; */
        }
        cout << endl;
    };
    return 0;
}

int IsingMatrix::initialize_spins(){
	Nspin = 0;
    for (int i = 0; i < NN; ++i){
        if (dis(gen) > 0.5){
            Nspin += 1;
            site[i].spin = true;
        } else {
            Nspin -= 1;
            site[i].spin = false;
        }
    };

    for (unsigned int i = 0; i < NN; ++i){
        site[i].n_bonds = 0;
    }

	Nbond = 0;
	for (unsigned int i = 0; i < NN; ++i){
		if (site[i].spin ^ site[i].right->spin)	{
			site[i].bond_right = false;
			site[i].n_bonds -= 1;
			site[i].right->n_bonds -= 1;
			Nbond -= 1;	
		} else {
			site[i].bond_right = true;
			site[i].n_bonds += 1;
			site[i].right->n_bonds += 1;
			Nbond += 1;
		}

		if (site[i].spin ^ site[i].down->spin)	{
			site[i].bond_down = false;
			site[i].n_bonds -= 1;
			site[i].down->n_bonds -= 1;
			Nbond -= 1;	
		} else {
			site[i].bond_down = true;
			site[i].n_bonds += 1;
			site[i].down->n_bonds += 1;
			Nbond += 1;
		}
	}		
    return 0;
}


IsingMatrix::IsingMatrix(int set_N, int set_Npick, int seed) 
  : N{set_N}, gen{std::mt19937(seed)}, dis {0.0, 1.0}, 
    Npick{set_Npick}, T_step_prior{-1.}, B_step_prior{-1.}
{ 
    /* cout << " seed: " << seed << " val: " << dis(gen) << endl; */
    NN = N*N;
    auto_correlation = new float[N/2 - 1];
    for (int i = 0; i < (N/2-1); ++i) auto_correlation[i] = 0;
	site = new lattice_site[NN];
	for (int i = 0; i < NN; ++i){
		indices.push_back(i);
	}

    n_spincorr = N/2 - (1 - N%2); // i.e. 4->1, 5->2, 6->2, 7->3, etc...
    spincorr = new float[n_spincorr];
    for (int i=0; i<n_spincorr; ++i) spincorr[i] = 0;

	// set up the links between the sites 
	//link left edge to right edge
	for (int i = 0; i < N; ++i){
		site[i*N].left        = &site[(i+1)*N-1];
		site[(i+1)*N-1].right = &site[i*N];
	}
	//link left to right, all other columns
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N - 1; ++j){
			site[i*N+j].right  = &site[i*N+j+1];
			site[i*N+j+1].left = &site[i*N+j];
		}
	}
	//link top edge to bottom edge
	for (int j = 0; j < N; ++j){
		site[j].up           = &site[N*(N-1) + j];
		site[N*(N-1)+j].down = &site[j];
	}
	//link bottom to top, all other rows
	for (int i = 0; i < N - 1; ++i){
		for (int j = 0; j < N; ++j){
			site[i*N+j].down = &site[(i+1)*N + j];
			site[(i+1)*N+j].up    = &site[i*N+j];
		}
	}

    /* for (int i{0}; i<NN;++i) site[i]->index = i; */

	initialize_spins(); 
}

IsingMatrix::~IsingMatrix(){ 
    delete[] site;
    delete[] auto_correlation;
    delete[] spincorr;
}

int IsingMatrix::calc_spincorr(){
    // calculate <ab>-<a><b> for all rows and columns
    
    // get row and col averages: <a> and <b>
    // 
    // note that <a> and <b> both equal the average spin of the matrix
    // (i.e. M)
    /* auto row_mean = new float[N]; */
    /* auto col_mean = new float[N]; */
    /* for (int i=0; i<N; ++i){ */
    /*     for (int j=0; j<N; ++j) { */
    /*         row_mean[i] += get_spin(i,j); */
    /*         col_mean[i] += get_spin(j,i); */
    /*     } */
    /*     row_mean[i] /= N; */
    /*     col_mean[i] /= N; */
    /* } */

    // reset the spin corr values
    for (int i=0; i<n_spincorr; ++i) spincorr[i] = 0;

    //calculate correlation of rows
    for (int row{0};row<N;++row){
        auto i_0 = row*N; // index of start of row
        for (int s{0};s<n_spincorr;++s){
            auto i_s = ((row+s+1)%N)*N; //index of start of row offset by s
            for (int col{0};col<N;++col){
                spincorr[s] += get_spin(i_0+col)*get_spin(i_s+col);
            }
        }
    }

    //calculate correlation of columns
    for (int col{0};col<N;++col){
        for (int s{0};s<n_spincorr;++s){
            int i_s { (col+s+1)%N }; //index of start of row offset by s
            for (int row{0};row<N;++row){
                int offset{row*N};
                spincorr[s] += get_spin(col+offset)*get_spin(i_s+offset);
            }
        }
    }
    
    const float M { get_M() };
    const float mean_sq { M*M };
    /* cout << " s: "; */
    for (int s{0};s<n_spincorr;++s)  spincorr[s] = spincorr[s]/(2.*NN) - mean_sq;
        /* if (s != 0) { */
        /*     cout << spincorr[s]<<"("<<(spincorr[s] == spincorr[s-1])<<") "; */
        /* } else { */
        /*     cout << spincorr[s] << " "; */
        /* } */
    /* cout << " s " << spincorr[1] << " " << spincorr[2] << " " << spincorr[2] - spincorr[1] << " : " << */ 
        /* (spincorr[1] == spincorr[2]) << endl; */

    return 0;
}


int IsingMatrix::calc_auto_correlation(){
    // get the row and column means
    vector<float> col_mean;
    vector<float> row_mean;
    for (int i = 0; i < N; ++i){
        vector<float> row_spin;
        vector<float> col_spin;
        for (int j = 0; j < N; ++j){
            row_spin.push_back( (float)get_spin(i,j) );
            col_spin.push_back( (float)get_spin(j,i) );
        }
        row_mean.push_back( accumulate(row_spin.begin(), row_spin.end(), 0.0)/N );
        col_mean.push_back( accumulate(col_spin.begin(), col_spin.end(), 0.0)/N );
    }
    /* cout << "row mean "; for (auto i : row_mean) cout << " " << i; cout << endl; */
    /* cout << "col mean "; for (auto i : col_mean) cout << " " << i; cout << endl; */
    for (int k = 1; k < N/2; ++k){
        float sum = 0;
        int n_entries = 0;
        for (int j = 1; j < N; ++j){
            for (int i = 0; i < N; ++i){
                int t = (j+k)%N;
                sum += (get_spin(j,i) - col_mean[i])*(get_spin(t,i)-col_mean[i]);
                sum += (get_spin(i,j) - row_mean[i])*(get_spin(i,t)-row_mean[i]);
                ++n_entries;
            }
        }
        sum /= (2*N*n_entries);
        /* sum /= (2*N); */
        auto_correlation[k-1] = sum;
    }
    return 0;
}

int IsingMatrix::nsteps(int n, float T, float B){
    for (int i = 0; i < n; ++i){
        step(T, B);
    }
    return 0;
}

void IsingMatrix::update_flip_prob(float T, float B){
    if ( (T!=T_step_prior) || (B!=B_step_prior)) {
        for (int spin : {0,1}) {
            for (int nbonds : {-4,-2,0,2,4}) {
                double deltaE { 2*(J*(nbonds)+B*(2.*spin-1)) };
                if (deltaE < 0) {
                    flip_prob[2+nbonds/2][spin] = 1;
                } else {
                    flip_prob[2+nbonds/2][spin] = exp(-1.0*deltaE/T);
                }
            }
        }
        T_step_prior = T;
        B_step_prior = B;
    }
}

int IsingMatrix::step(float T, float B_in){
    B = B_in;
    update_flip_prob(T,B);
    // determine the flip probabilities.
    
	// randomly pick Npick sites to potentially flip
	int n_flipped = 0;
	/* int num_random = Npick; */
	auto begin = indices.begin();
	auto end   = indices.end();
	size_t left = std::distance(begin, end);
	/* cout << " left " << left << endl; */
	/* cout << "begin" << endl; */
	/* cout << indices[0] << endl; */
	for (int i = 0; i < Npick; ++i){
	/* while(num_random--){ */
		auto r = begin;	
		float temp = (int)(dis(gen)*left);
		/* cout << " " << temp; */
		/* advance(r, (int)(dis(gen)*left)); // find a random number past begin(), */
		advance(r, temp); // find a random number past begin(),

		/* cout << " picked " << *r; */

        // the probability of swapping a site is calculated once in flip_prob
        // the value of -1 corresponds to guarantee that it will flip
		// the site r is allowed to flip if:
        //    1. deltaE < 0 or 2. exp(-1.0*deltaE/T) > rand()
		/* float deltaE = 2*(J*site[*r].n_bonds+B*(2*site[*r].spin-1)); */
		/* cout << " " << *r << " spin: " << site[*r].spin << " " << deltaE << endl; */
		/* cout << " deltaE: " << deltaE << " -> exp(-1*deltaE/T)=" << exp(-1.0*deltaE/T) << " "; */
        float prob = flip_prob[2+site[*r].n_bonds/2][site[*r].spin];
		if ( (prob == 1) || (dis(gen) < prob)) { //deltaE < 0.0 || dis(gen) < exp(-1.0*deltaE/T)){
			/* cout << " YES swap " << endl; */
			//move r to the front
			swap(*begin, *r);
			++begin;
			++n_flipped;
		} else {
			/* cout << " NO swap " << endl; */
		    // move r to the end, where it cannot get picked
			swap(*(end-1), *r);
		}
		--left;
		/* cout << " begin-end " << *begin << " " << *(begin + left -1) << endl; */
		/* cout << " v: "; for (auto k : indices) cout << " " << k; cout << endl; */
	}
	/* cout << endl; */
	/* cout << " indices: "; for (auto i : indices) cout << " " << i; cout << endl; */
	// now the first Npick numbers in indeces are
	// the sites that are allowed to flip.
	// check them all at once to see if they should flip,
	// and then flip them.
	/* cout << " index flipped: "; */
	for (int index = 0; index < n_flipped; ++index){	
		int i = indices[index];
		flip_site(site[i]);
	}
	return 0;	
}

// make a C-library
extern "C" // Tells the compiler to use C-linkage for this scope.
           // These are the functions available to ctypes in python
{
    /* void* c_new_IsingMatrix( void ) { return new(std::nothrow) IsingMatrix; } */
    void* newMatrix(int N, int n_flip, int seed_offset=0) { 
        return new IsingMatrix(N, n_flip, seed_offset);
    };
    int delMatrix( IsingMatrix* im) { delete im; return 0; }
    float get_E( IsingMatrix* im )  { return im->get_E(); }
    float get_M( IsingMatrix* im )  { return im->get_M(); }
    int   set_Npick( IsingMatrix* im, int val ) {
        return im->set_Npick( val ); }
    int   print_E( IsingMatrix* im )   { return im->print_E(); };
	int   print_Nspin( IsingMatrix* im ) { return im->print_Nspin(); };
	int   get_Nspin( IsingMatrix* im ) { return im->get_Nspin(); };
	int   get_Nbond( IsingMatrix* im ) { return im->get_Nbond(); };
	int   print_Nbond( IsingMatrix* im ) { return im->print_Nbond(); };
    int   get_spin( IsingMatrix* im, int i, int j) { return im->get_spin(i,j); };
    float get_site_flip_prob( IsingMatrix* im, int i, int j, float T, float B) 
          { return im->get_site_flip_prob(i,j,T,B); };
    int   set_spin( IsingMatrix* im, int i, int j, int k ) { return im->set_spin(i,j,k); };
    int   print_spins( IsingMatrix* im ) { return im->print_spins(); };     // defined
    int   print_bonds( IsingMatrix* im ) { return im->print_bonds(); };     // defined
    int   initialize_spins( IsingMatrix* im ){ return im->initialize_spins(); }; // defined
    int   step( IsingMatrix* im, float T, float B) { return im->step(T, B); };
    int   wolff( IsingMatrix* im, float T) { return im->wolff(T); };
    int   nsteps( IsingMatrix* im, int i, float T, float B) { return im->nsteps(i, T, B); };
    int   calc_auto_correlation( IsingMatrix* im ) { return im->calc_auto_correlation(); };
    int   calc_spincorr( IsingMatrix* im ) { return im->calc_spincorr(); };
    float auto_correlation( IsingMatrix* im, int i ){
        if (i >= im->get_N()/2+1){
            cout << "warning, there are only " 
                 << im->get_N()/2 
                 << " autocorrelation values available. " << endl;
            return -999.;
        }
        /* cout << "c i: " << i << " val: " << im->auto_correlation[i-1] << endl; */
        return im->auto_correlation[i-1]; 
    };
    float spincorr( IsingMatrix* im, int i) {
        if (i > im->n_spincorr) {
            cout << "n_spincorr " << im->n_spincorr << endl;
            cout << "warning, there are only " << im->n_spincorr
                 << " spin correlation values available." << endl;
            return -999.;
        }
        return im->spincorr[i];
    };
    int n_spincorr( IsingMatrix* im ) { return im->n_spincorr; };
};


int main(){
    cout << "starting main" << endl;
    IsingMatrix* z = reinterpret_cast<IsingMatrix*>(newMatrix(5,2));
    wolff(z, 2.5);
    cout << "*************" << endl;
    wolff(z, 2.5);
    return 0;
}
