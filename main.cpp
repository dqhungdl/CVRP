#include <bits/stdc++.h>
#include "ACO.hpp"
#include "BipartiteMatching.hpp"
using namespace std;

ofstream flog("log.txt");

struct location
{
	int64_t x,y,demand;
};

typedef pair<int64_t,int64_t> ii;
typedef pair<ii,ii> iii;

const int64_t inf=1e12;
const int64_t NumberofLoops=20;
const int64_t NumberofAttempts=3;
const int64_t MAXN=2000;
const int64_t Penalty=1e9;
const double EXP=0.2;

int64_t n,k,C;
int64_t LongestRoute;
int64_t Tmin=1,Tmax,T[MAXN+5][MAXN+5];
bool Used[MAXN+5];
vector<int64_t> rs[MAXN+5];
vector<ii> Route[MAXN+5];
vector<location> Location;
Problem instance;

// Rounding number
int64_t Round(double x)
{
	return int64_t(x+0.5);
}

// Get distance between (u1, v1) and (u2, v2)
int64_t Distance(int64_t u1,int64_t v1,int64_t u2,int64_t v2)
{
	return Round(sqrt((u1-u2)*(u1-u2)+(v1-v2)*(v1-v2)));
}

// Random in range [0; n/2]
int64_t NumberofWaitingActions(double n)
{
	double ratio=double(rand()%1001)/1000.0;
	return Round(n*ratio);
}

// Calculate penalty
int64_t SumPenalty(int64_t exceed)
{
	return exceed*(exceed+1)*Penalty;
}

// Calculate for reinforcement
int64_t SumT(int64_t id,int64_t u)
{
	int64_t sum=0;
	for(auto v:Route[id])
		sum+=T[u][v.first];
	return sum;
}

//******************************ACOTSP*******************************************************
// Optimize particular route
void ACOTSP()
{
	Solution solution;
	for(int64_t i=0;i<k;i++)
	{
		Tour tour;
		tour.id=i;
		for(auto tmp:Route[i])
		{
			tour.order.push_back(tmp.first);
			tour.tot_demand+=Location[tmp.first].demand;
		}
		solution.tours.push_back(tour);
	}
	ACO aco(&instance);
	aco.init_ACO();
	aco.fix_solution(solution);
	for(int64_t i=0;i<k;i++)
	{
		Route[i].clear();
		Tour TMP=solution.tours[i];
		for(auto tmp:TMP.order)
			Route[i].push_back(ii(tmp,0));
	}
}

void ACOTSP(int64_t i)
{
	vector<int> tour(Route[i].size());
	for(int64_t j=0;j<int64_t(Route[i].size());j++)
		tour[j]=Route[i][j].first;
	LS::ls_vector(tour,instance);
	for(int64_t j=0;j<int64_t(Route[i].size())-1;j++)
		Route[i][j].first=tour[j];
}
//*******************************************************************************************

//******************************Scoring particular solution**********************************
int64_t curminsumdist,minsumdist=1e18;
int64_t cost,overload;
bool Free[MAXN+5][MAXN+5];

void Statistic()
{
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			Free[i][j]=false;
	for(int i=0;i<k;i++)
		for(int j=0;j<int(Route[i].size());j++)
			for(int t=0;t<int(Route[i].size());t++)
				if(j!=t)
					Free[Route[i][j].first][Route[i][t].first]=true;
	int dup=0,dif=0;
	for(int i=0;i<k;i++)
		for(int j=0;j<int(rs[i].size());j++)
			for(int t=0;t<int(rs[i].size());t++)
				if(j!=t)
				{
					dup+=Free[rs[i][j]][rs[i][t]];
					dif+=1-Free[rs[i][j]][rs[i][t]];
				}
	cerr<<"Number of coincidences: "<<dup<<"\n";
	cerr<<"Number of differences: "<<dif<<"\n";
}

void PrintLog()
{
	for(int64_t i=0;i<70;i++)
		flog<<'-';
	flog<<"\n";
	for(int64_t i=0;i<k;i++)
	{
		flog<<"Route #"<<i+1<<": ";
		for(int64_t j=1;j<int64_t(rs[i].size())-1;j++)
			flog<<rs[i][j]<<' ';
		flog<<"\n";
	}
	flog<<"Cost "<<cost<<"\n";
}

int64_t Scoring()
{
	int64_t sumdist=0,_cost=0,_overload=0;
	for(int64_t i=0;i<k;i++)
	{
		int64_t sumcapacity=0;
		for(auto tmp:Route[i])
			sumcapacity+=Location[tmp.first].demand;
		if(sumcapacity>C)
		{
			sumdist+=SumPenalty(sumcapacity-C);
			_overload++;
		}
	}
	for(int64_t i=0;i<k;i++)
		for(int64_t j=1;j<int64_t(Route[i].size());j++)
		{
			int64_t u1=Location[Route[i][j-1].first].x,v1=Location[Route[i][j-1].first].y;
			int64_t u2=Location[Route[i][j].first].x,v2=Location[Route[i][j].first].y;
			sumdist+=Distance(u1,v1,u2,v2);
			_cost+=Distance(u1,v1,u2,v2);
		}
	if(minsumdist>sumdist)
	{
		Statistic();
		minsumdist=sumdist;
		cost=_cost;
		overload=_overload;
		for(int64_t i=0;i<k;i++)
		{
			rs[i].clear();
			for(auto tmp:Route[i])
				rs[i].push_back(tmp.first);
		}
		PrintLog();
	}
	return sumdist;
}
//*******************************************************************************************

//******************************Minimum Cost Bipartite Matching******************************
int64_t mm,nn,Start,Finish;
int64_t X[MAXN+5],Y[MAXN+5];								// Matching information
int64_t fx[MAXN+5],fy[MAXN+5];								// Kuhn-Munkres algorithm
int64_t Queue[MAXN+5],Trace[MAXN+5];						// BFS
int64_t c[MAXN+5][MAXN+5];									// Bipartite edges's weight
bool VisitedX[MAXN+5],VisitedY[MAXN+5];						// Vertices's mark
vector<int64_t> Side1,Side2,Capacity,ID;					// Two sides of bipartite graph

int64_t Get(int64_t i,int64_t j)
{
	return c[i][j]-fx[i]-fy[j];
}

void InitBipartiteGraph()
{
	Side1.clear();
	Side2.clear();
	Capacity.clear();
	ID.clear();
	// Init Side1
	for(int64_t i=0;i<2*k;i++)
	{
		Side1.push_back(Route[i].back().first);
		Capacity.push_back(Route[i].back().second);
		ID.push_back(i);
	}
	// Init Side2
	for(int64_t i=1;i<n;i++)
		if(Used[i]==false)
			Side2.push_back(i);						// Location i-th haven't been visited yet
	if(Side2.size()==0)
	{
		mm=0;
		return;
	}
	int64_t WaitingNumActions=NumberofWaitingActions(int64_t(Side1.size()));
	// Have no waiting action when starting the routes
	bool CheckStarted=false;
	for(int64_t i=0;i<int64_t(Side1.size());i++)
		if(Side1[i]!=0)
			CheckStarted=true;
	if(CheckStarted==false)
		WaitingNumActions=0;
	for(int64_t i=0;i<WaitingNumActions;i++)
		Side2.push_back(-1);						// Waiting action is numbered by -1
	while(Side1.size()>Side2.size())
		Side2.push_back(-1);						// For fully matching
	// Init bipartite graph's edges
	for(int64_t i=0;i<int64_t(Side1.size());i++)
		for(int64_t j=0;j<int64_t(Side2.size());j++)
		{
			if(Side2[j]==-1)
				c[i][j]=0;
			else
			{
				int64_t u1=Location[Side1[i]].x,v1=Location[Side1[i]].y;
				int64_t u2=Location[Side2[j]].x,v2=Location[Side2[j]].y;
				int64_t RequireCapacity=0;
				if(Side2[j]>=0)
					RequireCapacity=Location[Side2[j]].demand;
				c[i][j]=Distance(u1,v1,u2,v2)*SumT(ID[i],Side2[j]);//T[Side1[i]][Side2[j]];
				if(Capacity[i]+RequireCapacity>C)
					c[i][j]+=SumPenalty(Capacity[i]+RequireCapacity-C);
			}
		}
	// Init another information
	mm=Side1.size();
	nn=Side2.size();
	for(int64_t i=0;i<mm;i++)
	{
		X[i]=-1;
		fx[i]=0;
	}
	for(int64_t i=0;i<nn;i++)
	{
		Y[i]=-1;
		fy[i]=0;
	}
}

void Matching()
{
	for(int64_t i=0;i<nn;i++)
		Trace[i]=-1;
	Queue[1]=Start;
	int64_t l=1,r=1;
	while(l<=r)
	{
		int64_t i=Queue[l];
		l++;
		for(int64_t j=0;j<nn;j++)
			if(Trace[j]==-1&&Get(i,j)==0)
			{
				Trace[j]=i;
				if(Y[j]==-1)
				{
					Finish=j;
					return;
				}
				Queue[++r]=Y[j];
			}
	}
}

void SubX_AddY()
{
	for(int64_t i=0;i<mm;i++)
		VisitedX[i]=false;
	for(int64_t i=0;i<nn;i++)
		VisitedY[i]=false;
	VisitedX[Start]=true;
	for(int64_t i=0;i<nn;i++)
		if(Trace[i]!=-1)
		{
			VisitedX[Y[i]]=true;
			VisitedY[i]=true;
		}
	int64_t Delta=inf;
	for(int64_t i=0;i<mm;i++)
		if(VisitedX[i]==true)
			for(int64_t j=0;j<nn;j++)
				if(VisitedY[j]==false)
					Delta=min(Delta,Get(i,j));
	for(int64_t i=0;i<mm;i++)
		if(VisitedX[i]==true)
			fx[i]+=Delta;
	for(int64_t i=0;i<nn;i++)
		if(VisitedY[i]==true)
			fy[i]-=Delta;
}

void Enlarge()
{
	while(Finish!=-1)
	{
		int64_t x=Trace[Finish];
		int64_t next=X[x];
		X[x]=Finish;
		Y[Finish]=x;
		Finish=next;
	}
}

bool Hungarian()
{
	InitBipartiteGraph();
	if(mm==0)
		return false;
	// Hungarian algorithm
	for(int64_t i=0;i<mm;i++)
	{
		Start=i;
		Finish=-1;
		while(Finish==-1)
		{
			Matching();
			if(Finish==-1)
				SubX_AddY();
		}
		Enlarge();
	}
	// Modify routes
	for(int64_t i=0;i<mm;i++)
	{
		int64_t id=Side2[X[i]];
		if(id>=0)
		{
			Used[id]=true;
			Route[ID[i]].push_back(ii(id,Route[ID[i]].back().second+Location[id].demand));
		}
	}
	// Update longest route
	for(int64_t i=0;i<mm;i++)
		LongestRoute=max(LongestRoute,int64_t(Route[i].size()));
	return true;
}

bool cmp(vector<ii> x1,vector<ii> x2)
{
	return x1.back().second<x2.back().second;
}

void LastInitBipartiteGraph()
{
	sort(Route,Route+2*k,cmp);
	for(int64_t i=0;i<k;i++)
		for(int64_t j=k;j<2*k;j++)
		{
			int64_t u1=Location[Route[i].back().first].x,v1=Location[Route[i].back().first].y;
			int64_t u2=Location[Route[j].back().first].x,v2=Location[Route[j].back().first].y;
			c[i][j-k]=Distance(u1,v1,u2,v2);
			if(Route[i].back().second+Route[j].back().second>C)
				c[i][j-k]+=SumPenalty(Route[i].back().second+Route[j].back().second-C);
		}
	mm=nn=k;
	for(int64_t i=0;i<mm;i++)
	{
		X[i]=-1;
		fx[i]=0;
	}
	for(int64_t i=0;i<nn;i++)
	{
		Y[i]=-1;
		fy[i]=0;
	}
}

void LastHungarian()
{
	LastInitBipartiteGraph();
	// Hungarian algorithm
	for(int64_t i=0;i<mm;i++)
	{
		Start=i;
		Finish=-1;
		while(Finish==-1)
		{
			Matching();
			if(Finish==-1)
				SubX_AddY();
		}
		Enlarge();
	}
	// Matching routes
	for(int64_t i=0;i<mm;i++)
	{
		reverse(Route[X[i]+k].begin(),Route[X[i]+k].end());
		for(auto tmp:Route[X[i]+k])
			Route[i].push_back(tmp);
	}
}
//*******************************************************************************************

void ReadInput()
{
	string str;
	getline(cin,str);
	str+=" ";
	int64_t number=0;
	// Input n: Number of locatoins, k: Number of routes
	for(int64_t i=0;i<int64_t(str.size());i++)
		if('0'<=str[i]&&str[i]<='9')
			number=number*10+str[i]-'0';
		else
		if(number>0)
		{
			if(n==0)
				n=number;
			else
				k=number;
			number=0;
		}
	// Input capacity
	for(int64_t i=1;i<=4;i++)
		getline(cin,str);
	cin>>str>>str>>C;
	// Input the coordinate of all locations
	getline(cin,str);
	getline(cin,str);
	Location.resize(n);
	for(int64_t i=0;i<n;i++)
		cin>>number>>Location[i].x>>Location[i].y;
	// Input all demands
	getline(cin,str);
	getline(cin,str);
	for(int64_t i=0;i<n;i++)
		cin>>number>>Location[i].demand;
	// Read reinforcement learning data
	ifstream fi("ReinforcementLearning.txt");
	for(int64_t i=0;i<n;i++)
		for(int64_t j=0;j<n;j++)
			fi>>T[i][j];
	fi.close();
}

void InitInstance()
{
	instance.n=n;
	instance.distance=vector< vector<double> >(n,vector<double>(n,0));
	instance.nn_list=vector< vector<int> >(n,vector<int>(n-1,0));
	instance.points=vector<Point>(n);
	instance.demand.resize(n);
	for(int64_t i=0;i<n;i++)
	{
		instance.points[i].x=Location[i].x;
		instance.points[i].y=Location[i].y;
		instance.demand[i]=Location[i].demand;
	}
	// Compute distance between points
	for(int64_t i=0;i<n;i++)
		for(int64_t j=0;j<n;j++)
			instance.distance[i][j]=instance.points[i].round_Euclid_distance(instance.points[j]);
	// Compute nn_list
	for(int64_t i=0;i<n;i++) 
	{
		vector<pair<double,int> > nei;
		for(int64_t j=0;j<n;j++)
			if(j!=i)
				nei.push_back(make_pair(instance.distance[i][j],j));
		sort(nei.begin(),nei.end());
		for(int64_t j=0;j<n-1;j++)
			instance.nn_list[i][j]=nei[j].second;
	}
}

//*****************************Improvement***************************************************
bool ImprovedSwapLocations()
{
	bool CheckImproved=false;
	for(int64_t i=0;i<k;i++)
		for(int64_t j=i+1;j<k;j++)
			for(int64_t i1=1;i1<int64_t(Route[i].size())-1;i1++)
				for(int64_t j1=1;j1<int64_t(Route[j].size())-1;j1++)
				{
					swap(Route[i][i1].first,Route[j][j1].first);
					int64_t tmp=Scoring();
					if(curminsumdist<=tmp)
						swap(Route[i][i1].first,Route[j][j1].first);
					else
					{
						curminsumdist=tmp;
						vector<ii> Buffer1=Route[i],Buffer2=Route[j];
						ACOTSP(i);
						ACOTSP(j);
						tmp=Scoring();
						if(curminsumdist>tmp)
							curminsumdist=tmp;
						else
						{
							Route[i]=Buffer1;
							Route[j]=Buffer2;
						}
						CheckImproved=true;
					}
				}
	return CheckImproved;
}

bool ImprovedRoutes()
{
	bool CheckImproved=false;
	for(int64_t i=0;i<k;i++)
		for(int64_t j=i+1;j<k;j++)
			for(int64_t i1=0;i1<int64_t(Route[i].size())-1;i1++)
				for(int64_t j1=0;j1<int64_t(Route[j].size())-1;j1++)
				{
					vector<ii> NewRoute1,NewRoute2,Buffer1=Route[i],Buffer2=Route[j];
					for(int64_t t=0;t<=i1;t++)
						NewRoute1.push_back(Route[i][t]);
					for(int64_t t=j1+1;t<int64_t(Route[j].size());t++)
						NewRoute1.push_back(Route[j][t]);
					for(int64_t t=0;t<=j1;t++)
						NewRoute2.push_back(Route[j][t]);
					for(int64_t t=i1+1;t<int64_t(Route[i].size());t++)
						NewRoute2.push_back(Route[i][t]);
					Route[i]=NewRoute1;
					Route[j]=NewRoute2;
					ACOTSP(i);
					ACOTSP(j);
					int64_t tmp=Scoring();
					if(curminsumdist>tmp)
					{
						curminsumdist=tmp;
						CheckImproved=true;
					}
					else
					{
						Route[i]=Buffer1;
						Route[j]=Buffer2;
					}
					//debug(curminsumdist);
				}
	return CheckImproved;
}

void RandomACOTSP()
{
	int64_t nWhile=1000;
	while(nWhile--)
	{
		vector<ii> Buffer[MAXN+5];
		for(int64_t t=0;t<k;t++)
			Buffer[t]=Route[t];
		for(int64_t t=1;t<=3;t++)
		{
			int64_t i=rand()%k,j=rand()%k;
			int64_t i1=rand()%int64_t(Route[i].size()-2)+1;
			int64_t j1=rand()%int64_t(Route[j].size()-2)+1;
			if(i!=j)
				swap(Route[i][i1],Route[j][j1]);
		}
		ACOTSP();
		int64_t tmp=Scoring();
		if(curminsumdist>tmp)
		{
			cerr<<"Optimization at nWhile = "<<nWhile<<"\n";
			cerr<<"Current distance: "<<curminsumdist<<"\n";
			curminsumdist=tmp;
		}
		else
		for(int64_t t=0;t<k;t++)
			Route[t]=Buffer[t];
	}
}
//*******************************************************************************************
void Process()
{
	Tmax=n;
	// Init data
	LongestRoute=0;
	for(int64_t i=0;i<2*n;i++)
	{
		Route[i].clear();
		Route[i].push_back(ii(0,0));				// All routes start at depot(0)
	}
	for(int64_t i=0;i<n;i++)
		Used[i]=false;
	// Building k routes simultaneously
	while(Hungarian()==true){}
	LastHungarian();
	/*int nWhile=3;
	while(nWhile--)
	{
		curminsumdist=Scoring();
		while(ImprovedSwapLocations()==true||ImprovedRoutes()==true){}
		ACOTSP();
	}
	RandomACOTSP();*/
	Scoring();
}

bool Edge[MAXN+5][MAXN+5];

void ReinforcementLearning()
{
	for(int64_t i=0;i<n;i++)
		for(int64_t j=0;j<n;j++)
			Edge[i][j]=false;
	for(int64_t i=0;i<k;i++)
		for(int64_t j=0;j<int64_t(rs[i].size());j++)
			for(int64_t t=0;t<int64_t(rs[i].size());t++)
				if(j!=t)
					Edge[rs[i][j]][rs[i][t]]=true;
	for(int64_t i=0;i<n;i++)
		for(int64_t j=0;j<n;j++)
		{
			if(Edge[i][j]==true)
				T[i][j]-=Round(EXP*double(T[i][j]-Tmin));
			else
				T[i][j]-=Round(EXP*double(T[i][j]-Tmax));
			T[j][i]=T[i][j];
		}
	ofstream fo("ReinforcementLearning.txt");
	for(int64_t i=0;i<n;i++)
	{
		for(int64_t j=0;j<n;j++)
			fo<<T[i][j]<<' ';
		fo<<"\n";
	}
	fo.close();
}

void PrintOutput()
{
	for(int64_t i=0;i<k;i++)
	{
		cout<<"Route #"<<i+1<<": ";
		cout<<rs[i].size()<<' ';
		int64_t sumcapacity=0;
		for(auto tmp:rs[i])
		{
			cout<<tmp<<' ';
			if(tmp>0)
				sumcapacity+=Location[tmp].demand;
		}
		cout<<"\n"<<"Capacity #"<<i+1<<": "<<sumcapacity;
		if(sumcapacity>C)
			cout<<"  Overload!";
		cout<<"\n";
	}
	cout<<"Cost: "<<cost<<"\n";
	cout<<"Number of overload vehicles: "<<overload;
}

void PrintOutput2()
{
	ofstream fo("output.txt");
	for(int64_t i=0;i<k;i++)
	{
		fo<<"Route #"<<i+1<<": ";
		for(int64_t j=1;j<int64_t(rs[i].size())-1;j++)
			fo<<rs[i][j]<<' ';
		fo<<"\n";
	}
	fo<<"Cost "<<cost;
	fo.close();
}

int main()
{
	srand(time(0));
	freopen("TEST.txt","r",stdin);
	freopen("TEST.ans","w",stdout);
	for(int64_t T=1;T<=NumberofLoops;T++)
	{
		cerr<<NumberofLoops-T+1<<" loop(s) remaining\n";
		ReadInput();
		InitInstance();
		for(int64_t i=1;i<=NumberofAttempts;i++)
		{
			Process();
			cerr<<"\tAttemp #"<<i<<":  Cost: "<<cost<<",  Overload: "<<overload<<"\n";
		}
		cerr<<"Current min length: "<<cost<<"\n";
		cerr<<"Current min overload: "<<overload<<"\n";
		//ReinforcementLearning();
	}
	PrintOutput();
	PrintOutput2();
	PAUSE();
}