#include <bits/stdc++.h>
using namespace std;

struct location
{
	int64_t x,y,demand;
};

typedef pair<int64_t,int64_t> ii;

const int64_t inf=1e12;
const int64_t NumberofAttempts=50;
const int64_t MAXN=2000;
const int64_t Penalty=1e9;
const double EXP=0.5;

int64_t n,k,C;
int64_t depotX,depotY;
int64_t Tmin=1,Tmax,T[MAXN+5][MAXN+5];
int64_t NumberofFines=0;
bool Used[MAXN+5];
vector<int64_t> rs[MAXN+5];
vector<ii> Route[MAXN+5];
vector<location> Location;
multiset<int64_t> setDemands;

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

// Random in range [n/2.5; n/1.5]
int64_t NumberofWaitingActions(double n)
{
	double ratio=double(rand()%101+150)/100.0;
	return Round(n/ratio);
}

// Scoring particular solution
ii minsumdist=ii(1e12,0);

void Scoring()
{
	ii sumdist=ii(0,0);
	for(int64_t i=0;i<k;i++)
	{
		int64_t sumcapacity=0;
		for(auto tmp:Route[i])
			if(tmp.first>0)
				sumcapacity+=Location[tmp.first-1].demand;
		if(sumcapacity>C)
			sumdist.second++;
	}
	for(int64_t i=0;i<k;i++)
	{
		for(int64_t j=1;j<int64_t(Route[i].size());j++)
		{
			int64_t u1,v1,u2,v2,id1=Route[i][j-1].first-1,id2=Route[i][j].first-1;
			if(id1==-1)
			{
				u1=depotX;
				v1=depotY;
			}
			else
			{
				u1=Location[id1].x;
				v1=Location[id1].y;
			}
			if(id2==-1)
			{
				u2=depotX;
				v2=depotY;
			}
			else
			{
				u2=Location[id2].x;
				v2=Location[id2].y;
			}
			sumdist.first+=Distance(u1,v1,u2,v2);
		}
	}
	if(minsumdist.first+Penalty*minsumdist.second>sumdist.first+Penalty*sumdist.second)
	{
		minsumdist=sumdist;
		for(int64_t i=0;i<k;i++)
		{
			rs[i].clear();
			for(auto tmp:Route[i])
				rs[i].push_back(tmp.first);
		}
	}
}

//******************************Minimum Cost Bipartite Matching******************************
int64_t mm,nn,Start,Finish;
int64_t X[MAXN+5],Y[MAXN+5];								// Matching information
int64_t fx[MAXN+5],fy[MAXN+5];								// Kuhn-Munkres algorithm
int64_t Queue[MAXN+5],Trace[MAXN+5];						// BFS
int64_t c[MAXN+5][MAXN+5];									// Bipartite edges's weight
bool VisitedX[MAXN+5],VisitedY[MAXN+5];						// Vertices's mark
vector<int64_t> Side1,Side2,RemainCapacity,ID;				// Two sides of bipartite graph

int64_t Get(int64_t i,int64_t j)
{
	return c[i][j]-fx[i]-fy[j];
}

void InitBipartiteGraph()
{
	Side1.clear();
	Side2.clear();
	RemainCapacity.clear();
	ID.clear();
	// Init Side1
	for(int64_t i=0;i<k;i++)
		if(Route[i].size()==1)
		{
			Side1.push_back(-1);					// Depot is numbered by -1
			RemainCapacity.push_back(C);			// Init capacity is C
			ID.push_back(i);						// ID of route
		}
		else
		if(Route[i].back().first!=0)
		{
			Side1.push_back(Route[i].back().first-1);
			RemainCapacity.push_back(Route[i].back().second);
			ID.push_back(i);
		}
	// Init Side2
	for(int64_t i=0;i<n;i++)
		if(Used[i]==false)
			Side2.push_back(i);						// Location i-th haven't been visited yet
	while(Side1.size()>Side2.size())
		Side2.push_back(-1);						// Go back to depot -1
	int64_t WaitingNumActions=NumberofWaitingActions(int64_t(Side1.size()));
	for(int64_t i=0;i<WaitingNumActions;i++)
		Side2.push_back(-2);						//Waiting action is numbered by -2
	// Init bipartite graph's edges
	for(int64_t i=0;i<int64_t(Side1.size());i++)
		for(int64_t j=0;j<int64_t(Side2.size());j++)
		{
			if(Side2[j]==-2)
				c[i][j]=0;
			else
			{
				int64_t u1,v1,u2,v2;
				if(Side1[i]==-1)
				{
					u1=depotX;
					v1=depotY;
				}
				else
				{
					u1=Location[Side1[i]].x;
					v1=Location[Side1[i]].y;
				}
				if(Side2[j]==-1)
				{
					u2=depotX;
					v2=depotY;
				}
				else
				{
					u2=Location[Side2[j]].x;
					v2=Location[Side2[j]].y;
				}
				int64_t RequireCapacity=0;
				if(Side2[j]>=0)
					RequireCapacity=Location[Side2[j]].demand;
				c[i][j]=Distance(u1,v1,u2,v2)*T[Side1[i]+1][Side2[j]+1];
				if(RemainCapacity[i]>=RequireCapacity&&RemainCapacity[i]>=*setDemands.begin()&&Side2[j]==-1)
					c[i][j]+=Penalty;
				else
				if(RemainCapacity[i]<RequireCapacity)
					c[i][j]+=inf;
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
		if(id==-1)
			Route[ID[i]].push_back(ii(0,0));
		else
		if(id>=0)
		{
			Used[id]=true;
			setDemands.erase(setDemands.find(Location[id].demand));
			if(Route[ID[i]].back().second<Location[id].demand)
				Route[ID[i]].push_back(ii(id+1,0));
			else
				Route[ID[i]].push_back(ii(id+1,max(Route[ID[i]].back().second-Location[id].demand,0LL)));
		}
	}
	return true;
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
	int64_t sumdemand=0;
	for(int64_t i=0;i<n;i++)
	{
		cin>>number>>Location[i].demand;
		sumdemand+=Location[i].demand;
	}
	// Input the depot's coordinate
	getline(cin,str);
	getline(cin,str);
	cin>>depotX>>depotY;
	// Read reinforcement learning data
	ifstream fi("ReinforcementLearning.txt");
	for(int64_t i=0;i<=n;i++)
		for(int64_t j=0;j<=n;j++)
			fi>>T[i][j];
	fi.close();
}

void Process()
{
	// Init data
	for(int64_t i=0;i<n;i++)
	{
		setDemands.insert(Location[i].demand);
		Route[i].clear();
		Route[i].push_back(ii(0,C));				// All routes start at depot(0)
	}
	for(int64_t i=0;i<n;i++)
		Used[i]=false;
	// Building k routes simultaneously
	while(Hungarian()==true){}
	Scoring();
}

bool Edge[MAXN+5][MAXN+5];

void ReinforcementLearning()
{
	Tmax=n;
	for(int64_t i=0;i<=n;i++)
		for(int64_t j=0;j<=n;j++)
			Edge[i][j]=false;
	for(int64_t i=0;i<k;i++)
		for(int64_t j=1;j<int64_t(rs[i].size());j++)
			Edge[rs[i][j-1]][rs[i][j]]=true;
	for(int64_t i=0;i<=n;i++)
		for(int64_t j=0;j<=n;j++)
			if(Edge[i][j]==true)
				T[i][j]-=Round(EXP*double(T[i][j]-Tmin));
			else
				T[i][j]-=Round(EXP*double(T[i][j]-Tmax));
	ofstream fo("ReinforcementLearning.txt");
	for(int64_t i=0;i<=n;i++)
	{
		for(int64_t j=0;j<=n;j++)
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
				sumcapacity+=Location[tmp-1].demand;
		}
		cout<<"\n"<<"Capacity #"<<i<<": "<<sumcapacity;
		if(sumcapacity>C)
			cout<<"  Overload!";
		cout<<"\n";
	}
	cout<<"Cost: "<<minsumdist.first<<"\n";
	cout<<"Number of overload vehicles: "<<minsumdist.second;
}

int main()
{
	srand(time(0));
	freopen("TEST.txt","r",stdin);
	freopen("TEST.ans","w",stdout);
	for(int64_t i=1;i<=50;i++)
	{
		ReadInput();
		for(int64_t i=1;i<=NumberofAttempts;i++)
		{
			debug(i);
			Process();
		}
		ReinforcementLearning();
	}
	PrintOutput();
}