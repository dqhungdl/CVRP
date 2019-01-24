#include <bits/stdc++.h>
using namespace std;

ifstream fi("output.txt");

struct location
{
	int x,y,demand;
};

const int MAXN=2000;

int n,k,C;
vector<location> Location;

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
}

int Distance(int u1,int v1,int u2,int v2)
{
	double tmp=sqrt((u1-u2)*(u1-u2)+(v1-v2)*(v1-v2));
	return int(tmp+0.5);
}

void ReadOutput()
{
	int num,id,sumdist=0;
	string str;
	for(int i=0;i<k;i++)
	{
		fi>>str>>str>>num;
		vector<int> Route;
		while(num--)
		{
			fi>>id;
			Route.push_back(id);
		}
		int sumdemand=0;
		for(int i=1;i<int(Route.size());i++)
		{
			sumdemand+=Location[Route[i]].demand;
			sumdist+=Distance(Location[Route[i-1]].x,Location[Route[i-1]].y,Location[Route[i]].x,Location[Route[i]].y);
		}
		cout<<"Capacity #"<<id+1<<": "<<sumdemand<<"\n";
		fi>>str>>str>>str;
	}
	cout<<"Cost "<<sumdist;
}

int main()
{
	freopen("TEST.txt","r",stdin);
	ReadInput();
	ReadOutput();
	PAUSE();
}