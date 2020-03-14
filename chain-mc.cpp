#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib> 
#include <ctime>

using namespace std;
typedef std::numeric_limits< double > dbl;

int mod(int a, int b) //a proper mod function //FIX THIS!
{
    int ret = a%b;
    if(ret<0)
        ret += b;
    return ret;
}
double urand(void) //random # in [0,1)
{
	return ((double)rand())/((double)RAND_MAX);
}
class Lattice
{
public:
	int Nx, Ny, N;
	int Nb;
	int *state;
	int **bSites;
	Lattice() //Empty lattice. 
	{
		N = 0;
	}
	Lattice(int N) // length N chain
	{
		this->N = N;
		initState();
		initBonds();
	}
	Lattice(int N, int numUp, bool what) //lattice with numUp spins up
	{
		//boolean does nothing but distinguish the constructor. it's silly
		this->N = N;
		initState(numUp);
		initBonds();
	}
	Lattice(int Nx, int Ny) // Nx by Ny lattice
	{
		this->Nx = Nx; this-> Ny = Ny;
		this->N = Nx*Ny;
		initState();
		init2DBonds();

	}
	~Lattice(void)
	{
		//cout << "Deconstruct Lattice" << endl;
		for(int i=0;i<Nb;i++)
			delete [] bSites[i];
		delete [] bSites;
		delete [] state;
	}
	void printBonds(void)
	{
		for(int b=0;b < Nb; b++)
			cout << bSites[b][0] << "," << bSites[b][1] << endl;
	}
private:
	void initState(void) //generate random initial state
	{
		state = new int[N];
		int sum;
		do{
			sum = 0;
			for(int i=0;i<N;i++)
			{
				state[i] = (rand()%2)*2-1;
				sum += state[i];
			}
		}while(abs(sum)==N); //exclude all up/down configurations
	}
	void initState(int numUp) //state with numUp spins up (not random)
	{
		state = new int[N];
		for(int i=0;i<N;i++)
		{
			if(i<numUp)
				state[i] = 1;
			else
				state[i] = -1;
		}
	}
	void initBonds(void) //define bonds linking sites on chain
		{
		Nb = N;
		bSites = new int*[Nb];
		for(int i=0; i<Nb; i++)
			{
				bSites[i] = new int[2];
				bSites[i][0] = i;
				bSites[i][1] = (i+1)%Nb;
			}
	
		}
	void init2DBonds(void) //bonds on 2D lattice
	{
		Nb = 2*N;
		bSites = new int*[Nb];
		for(int i=0; i<N; i++)
		{
			bSites[i] = new int[2];
			bSites[i+N] = new int[2];
			bSites[i][0] = i;
			bSites[i][1] = (i+1)%Nx + (i/Nx)*Nx;
			bSites[i+N][0] = i;
			bSites[i+N][1] = (i+Nx)%N;
		}
	
	}
};

class OpList
{
private:
	

public:
	int nMax;
	int M; //max word length
	int n; //word length
	int* word; //word stores operators and their location on lattice by -1 for identity, and (bond #)*opType for others
	int* links;
	enum opType {iden=-1,diag=0,oDiag=1};
	static const int opTypes = 2;
	char names[2];
	
	OpList(int M)
	{
		this->M = M;
		word = new int[M];
		for(int p=0;p<M;p++)
			word[p]=-1; //set all operators to identity
		n = 0;
		nMax = n;
		links = new int[M*4];
		names[0] = 'D'; names[1] = 'O';
	}
	~OpList()
	{
		//cout << "Deconstruct opList" << endl;
		//delete [] word;
		//delete [] links;
	}
	void adjustLength(void) //make M larger if n is too big
	{
		if(float(nMax) > 3*float(M)/4)
		{
			int oldM = M;
			M = int(float(n)*4/3+1);
			int* newWord = new int[M];
			for(int p=0;p<M;p++)
			{
				if(p<oldM)
					newWord[p] = word[p];
				else
					newWord[p] = iden;
			}
			delete [] word;
			word = newWord;
		}
	}
	bool isType(int letter, opType op) //check if letter in word is operator op
	{
		if((letter == op) || (op == iden))
			return (letter == (int)op);
		return (letter%opTypes == (int)(op));
	}
	int getBond(int letter) //returns bond the operator is on at letter
	{
		return letter/2;
	}
	void showWord(Lattice &latt) //show word + lattice changes
	{
		int* curState;
		curState = new int[latt.N];
		for(int i=0;i<latt.N;i++)
			curState[i] = latt.state[i];
		//std::copy(std::begin(latt.state),std::end(latt.state),std::begin(curState));
		int N = latt.N; int Nb = latt.Nb;
		int** bSites = latt.bSites;
		string dots = "";
		for(int i=0;i<N;i++)
		{
			if(curState[i]>0)
				dots += "x ";
			else
				dots += ". ";
		}
		cout << dots << endl;
		for(int p=0; p<M; p++)
		{
			for(int n=0;n<Nb;n++)
				if(isType(word[p],iden))
					cout << "  ";
				else
					if(n == word[p]/opTypes)
						cout << " " << names[word[p]%opTypes];
					else
						cout << "  ";
			if(isType(word[p],oDiag))
			{
				
				int bond = word[p]/opTypes;
				curState[bSites[bond][0]] *= -1;
				curState[bSites[bond][1]] *= -1;
			}
			dots = "";
			for(int i=0;i<N;i++)
			{
				if(curState[i]>0)
					dots += "x ";
				else
					dots += ". ";
			}
			cout << endl << dots << endl;
			//delete [] curState;
		}
	}
	void constructLinks(Lattice &latt) //construct linked list of operators
	{
		int** bSites = latt.bSites; //for convenient reference
		int len = M*4; //4 links per operator
		links = new int[len];
		for(int v = 0;v<len;v++) //init all links as unlinked
			links[v] = -1;
		for(int p0 = 0;p0<M;p0++) //loop thru operator string
			if(!isType(word[p0],iden)) //ignore identity
				for(int l=0;l<4;l++) //loop thru links at current operator
				{
					int v0 = 4*p0 + l; //current position in linked list
					if(links[v0]== -1) //if current link hasn't been linked, find the operator linked to it
					{
						int p = p0;
						int dp = (l/2)*2-1; // -1 for l=0,1 & +1 for l=3,4 (moving up or down list)
						bool stop = false;
						int l1, v1;
						while(!stop) //move along list until we find a friend for our operator (which can be itself!)
						{
							p = mod(p+dp,M); //move up/down list
							if(!isType(word[p],iden)) //if operator here is not identity
							{
								int b0 = getBond(word[p0]); //bond of our first operator
								int b1 = getBond(word[p]);  //bond of new potential friend
								if(bSites[b0][l%2]==bSites[b1][0])//connect on left
								{
									l1 = 2-2*(l/2); //corresponding link on friend
									stop = true; //stop after finding friend
								}	
								if(bSites[b0][l%2]==bSites[b1][1]) //on right
								{
									l1 = 3-2*(l/2);
									stop = true;
								}
							}
						}
						v1 = 4*p + l1; //position in linked list for friend
						//link operator to its friend and vice versa
						links[v0] = v1;
						links[v1] = v0;
					}
				}
	}
	void printLinks(void)
	{
		cout << "v: ";
		for(int v=0;v < 4*M;v++)
			cout << links[v] << " ";
		cout << endl;
	}
};
void checkPeriodicity(int* curState, Lattice &latt) //check final state = initial state (for debug)
{
	for(int i=0; i<latt.N; i++)
		if(curState[i]!=latt.state[i])
			cout << "BAD NEWS BEARS EVERYONE" << endl;
}
OpList diagonalStep(Lattice &latt, OpList &opList,double beta)
{
	int Nb = latt.Nb; int M = opList.M; int n = opList.n;
	int* curState;
	curState = new int[latt.N];
	for(int i=0;i<latt.N;i++)
		curState[i] = latt.state[i]; //current state is copied so we can evolve it to get state at each euclidean time step
	int** bSites = latt.bSites;
	int* word = opList.word;
	//probabilities for operator insertion/deletion
	double pDiag = min(Nb*beta/2/(double(M-n)),(double)1);
	double pId = min(2*((double)(M-n+1))/(Nb*beta),(double)1);
	//loop through operator list
	for(int p=0;p<M;p++)
	{
		if(opList.isType(word[p],opList.iden)) //if identity
		{
			int bond = rand()%Nb; //random bond
			//accept bond only if it has antiparallel spins
			if(curState[bSites[bond][0]]*curState[bSites[bond][1]]<0)
				if(urand() < pDiag)
				{
					word[p] = opList.opTypes*bond; //change operator at p to diagonal at location bond
					n++ ; //increment count of non identity operators
				}
		}
		else
		{
			if(opList.isType(word[p],opList.diag)) //if diagonal
			{
				if(urand() < pId)
				{
					word[p] = opList.iden;//change operator at p to identity
					n-- ; //decrement count of non identity operators
				}
			
			}
			else //if operator neither identity nor diagonal (ie off diagonal)
			{
				int bond = word[p]/opList.opTypes; //i don't understand how this gives you the right bond
				//flip spins on off diagonal bond
				curState[bSites[bond][0]] *= -1;
				curState[bSites[bond][1]] *= -1;
			}
		}
	}
	//checkPeriodicity(curState,latt);
	opList.n = n; //update number of operators
	if (n>opList.nMax) //increase max number of operators if n exceeds it
		opList.nMax = n;
	return opList; //return updated operator list
}
void doLoop(int vi,Lattice &latt, OpList &opList) //flip loops, mark operators for change
{
	int* links = opList.links;
	int* word = opList.word;
	int* state = latt.state; //why is this by reference? 
	int** bSites = latt.bSites;
	int update = -1-(rand()%2); //-2 to update, -1 not to update loop
	int v0 = vi;
	int v2 = vi^1;
	do{
		int v1 = links[v0]; //follow link to next bond
		v2 = v1^1; //xor to get other site on bond
		links[v0] = update; //overwrite link values with flip choice
		links[v1] = update;
		if(update == -2) // flip spins in loop if update
		{
			int vhigh = max(v0,v1);
			int vlow = min(v0,v1);
			//flip spins in initial state if loop crosses it
			if(vhigh%4 > vlow%4) //if loop crosses initial state (loops around)
			{
				int bond = word[v0/4]/opList.opTypes;
				state[bSites[bond][v0%2]] *= -1; //get index of spin, flip spins in initial state
			}
		}
		v0 = v2;
	}while(v2 != vi);
}
void loopStep(Lattice &latt, OpList &opList)
{
	opList.constructLinks(latt);
	int* links = opList.links;
	int* word = opList.word;
	int M = opList.M;
	//loop thru linked list and create loops
	for(int v=0;v<4*M;v++) 
		if(links[v] > -1) //start loops on unvisted links
			doLoop(v,latt,opList);
	for(int p=0;p<M;p++)
		if(!opList.isType(word[p],opList.iden)) //skip identity operators
		{
			int p4 = p*4;
			int flip = links[p4] + links[p4+1] + links[p4+2] + links[p4+3];
			if(flip == -6) // nonflip on one side, flip on the other
				word[p] = word[p]^1; //change diagonal to off diagonal and vice versa
		}
	delete [] opList.links;
}
double expectationE(Lattice &latt,double T,int itSteps,double conv)
{
	double beta = 1/T;
	int sumN = 0;
	int curStep = 0;
	OpList opList(latt.N);
	bool startAvg = false;
	double checkAvg = -100;
	double curAvg = 0;
	bool stop = false;
	for(int it=0; it<itSteps && !stop; it++)
	{
		opList = diagonalStep(latt,opList,beta);
		//opList.showWord(latt);
		//cout << "-------------" << endl;
		loopStep(latt,opList);
		//opList.showWord(latt);
		opList.adjustLength();
		
		if(it>itSteps/4)
		{
			curStep++ ;
			sumN += opList.n;
			curAvg = ((double)(sumN))/((double)curStep);
			if(it%1000==0)
				if(abs(curAvg-checkAvg)<conv)
					stop = true;
				else
					checkAvg = curAvg;
		}
		
		if(it%(itSteps/10)==0)
		{
			cout << "M= " << opList.M << ", n= " << opList.n;
			cout << ", avg= " << curAvg << endl;
		}
	}
	cout << curStep << " steps" << endl;
	return curAvg/beta;
}
double singleTemp(int dim,int N, double T, int itSteps, int states, double prec)
{
	double avgE = 0;
	for(int i=0; i<states; i++)
	{
		Lattice* latt;
		//delete latt;
		if(dim == 2)
			latt = new Lattice(N,N);
		else
			latt = new Lattice(N);
		avgE += expectationE(*latt,T,itSteps,prec);
		delete latt;
	}
	avgE = avgE/((double)states);
	return avgE;
}
bool readInput(char* fileName, double* param,
	int numParam, int lineSkip)
{
	ifstream inFile(fileName);
    if (!inFile.is_open())
    {
        cout << "Unable to open " << fileName << endl;
        return false;
    }
    for(int i = 0;i < lineSkip && !inFile.eof();i++)
    {
        string line;
        getline(inFile,line);
    }
	int i;
    for(i = 0;i<numParam && !inFile.eof();i++)
    {
        inFile >> param[i];
    }
    //abort if not all values are found, or file ends
    if(i<numParam || inFile.eof())
    {
        cout << "Improper input" << endl;
        return false;
    }
	return true;
}
int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Enter input file" << endl;
		cout << "Format: dim Steps States Precision nSizes nTemps" << endl;
		return 0;
	}
	const int nParam = 6;
	int nSizes;
	int nTemps;
	char* inFileName = argv[1];
	char* outFileName;
	ofstream outFile;
	if(argc > 2)
	{
		outFileName = argv[2];
		outFile.open(outFileName);
		if(!outFile.is_open())
		{
			cout << "Cannot open " << outFileName << endl;
			return 0;
		}
		outFile.close();
	}
	else
	{
		cout << "No output file" << endl;
		return 0;
	}
	cout.precision(dbl::digits10); //display double digits
	srand(time(0)); //new random seed
	//read in parameters
	double param[nParam];
	bool readFile = readInput(inFileName,param,nParam,0);
	if(!readFile) return 0;
	int dim = param[0]; int itSteps = param[1]; int states = param[2];
	double prec = param[3]; nSizes = param[4]; nTemps = param[5];
	//read in lattice sizes
	double *lattSizes;
	lattSizes = new double[nSizes];
	readFile = readInput(inFileName,lattSizes,nSizes,1);
	if(!readFile) return 0;
	//read in temperatures
	double *temps;
	temps = new double[nTemps];
	readFile = readInput(inFileName,temps,nTemps,2);
	if(!readFile) return 0;
	//initialize output array
	double *energy;
	energy = new double[nTemps*nSizes];
	cout << "dim = " << dim << ", ";
	cout << itSteps << " steps with precision " << prec << endl;
	cout << "nSizes = " << nSizes << ", nTemps = " << nTemps << endl;
	for(int i=0;i<nSizes;i++)
	{
		cout << "N = " << (int)lattSizes[i] << endl;
		for(int j=0;j<nTemps;j++)
		{
			cout << "T = " << temps[j] << endl;
			double E = singleTemp(dim,(int)lattSizes[i],temps[j],
									itSteps, states, prec);
			if(nSizes > nTemps)
				energy[i+j*nSizes] = E;
			else
				energy[j+i*nTemps] = E;
		}
		
	}
	cout << "argc: " <<  argc << endl;
	if(argc > 2)
	{
		cout << "Writing to " << outFileName << endl;
		outFile.open(outFileName);
		outFile.setf(ios::fixed,ios::floatfield);
		outFile.precision(dbl::digits10);
		for(int i=0;i<min(nTemps,nSizes);i++)
		{
			for(int j=0;j<max(nTemps,nSizes);j++)
			{
				outFile << energy[j+i*((int)max(nSizes,nTemps))] << "	";
			}
			outFile << endl;
		}
		outFile.close();
	}
	cout << "dim = " << dim << ", ";
	cout << itSteps << " steps with precision " << prec << endl;
	cout << "nSizes = " << nSizes << ", nTemps = " << nTemps << endl;
	delete [] lattSizes;
	delete [] temps;
	return 0;
}
