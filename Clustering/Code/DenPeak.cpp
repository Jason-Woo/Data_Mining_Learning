/*
datasets\\Aggregation_cluster=7.txt
datasets\\flame_cluster=2.txt
datasets\\Jain_cluster=2.txt
datasets\\Pathbased_cluster=3.txt
datasets\\Spiral_cluster=3.txt
*/
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<cmath>
#include<sstream>
#include<algorithm>

using namespace std;

struct point 
{
    double x;
    double y;
    int label;
    int id;
    bool center;
};

struct node
{
    int id;
    double rho;
    double delta;
};

class Cluster
{
    public:
        Cluster(string filename);
        int fclust();
        void getdist();
        void getdc();
        void getrho();
        void getdelta();
        void assign();
        void Output(string filename);
    private:
        vector<point> data;
        double dc;                         
        vector<vector<double> > t_dist; 
        double t_neighbor_rate;
        double maxdist;

        vector<double> t_rho;
        vector<node> t_orderrho;
        vector<double> t_delta;

        vector<int> t_neighbor;
        double t_maxrho;
        double t_minrho;

        int t_clusterNum;
        vector<int> cl;
        vector<int> icenter_of_class;    
};


Cluster::Cluster(string filename)
{
    t_neighbor_rate = 0.07;
    ifstream fin(filename);
    if(! fin)
    {
        cout << "OPEN FILE ERROR" << endl;
        return;
    }

    string buffer;
    int input_cnt = 0;
    while(getline(fin, buffer))
    {
        point p;
		stringstream ss(buffer);
        string tmp;
        ss >> tmp;
        double tmp_num = stod(tmp);
        p.x = tmp_num;
        ss >> tmp;
        tmp_num = stod(tmp);
        p.y = tmp_num;
        p.label = -1;
        p.id = input_cnt;
        p.center = false;
        data.push_back(p);
        input_cnt ++;
    }

    t_dist.resize(data.size());
    for(int i = 0 ; i < data.size(); i ++) t_dist[i].resize(data.size());
    t_rho.resize(data.size());
    t_orderrho.resize(data.size());
}

void Cluster::getdist() 
{
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = i; j < data.size(); j++)
		{
			double distij = sqrt(pow(data[i].x - data[j].x, 2) + pow(data[i].y - data[j].y, 2));
			t_dist[i][j] = distij;
			t_dist[j][i] = distij;
		}
	}
}

void Cluster::getdc()
{
	vector<double> alldist;
	for(int i = 0; i < t_dist.size() - 1; i ++)
	{
		for(int j = i + 1; j < t_dist.size(); j ++)
		{
			alldist.push_back(t_dist[i][j]);
		}
	}
	sort(alldist.begin(), alldist.end());
	dc = alldist[static_cast<int>(alldist.size()) * t_neighbor_rate];
	maxdist = alldist[alldist.size() - 1];
}

bool comp(node x, node y)
{
	return x.rho > y.rho;
}

void Cluster::getrho()
{	
	for(int i = 0; i < t_dist.size() - 1; i ++)
	{
		for(int j = i + 1; j < t_dist.size(); j ++)
		{
			double distij = t_dist[i][j];
			t_rho[i] = t_rho[i] + exp(-1 * (pow(t_dist[i][j] / dc, 2)));
			t_rho[j] = t_rho[j] + exp(-1 * (pow(t_dist[i][j] / dc, 2)));
		}
	}
	
	for(int i = 0; i < t_rho.size(); i ++)
	{
		t_orderrho[i].id = i;
		t_orderrho[i].rho = t_rho[i];
	}
	sort(t_orderrho.begin(), t_orderrho.end(), comp);
}

void Cluster::getdelta()
{
	t_delta.resize(data.size(), 0);
	t_neighbor.resize(data.size(), -1);
	t_delta[t_orderrho[0].id] = -1.0;
	t_neighbor[t_orderrho[0].id] = -1;
	
	for (int i = 0; i < data.size(); i ++)
	{
		t_delta[t_orderrho[i].id] = maxdist;
		for(int j = 0; j < i; j ++)
		{
			if(t_dist[t_orderrho[i].id][t_orderrho[j].id] < t_delta[t_orderrho[i].id])
			{
				t_delta[t_orderrho[i].id] = t_dist[t_orderrho[i].id][t_orderrho[j].id];
				t_neighbor[t_orderrho[i].id] = t_orderrho[j].id;
			}
		}
	}
	t_delta[t_orderrho[0].id] = maxdist;
}

void Cluster::assign()
{
	double t_rhoTH = t_maxrho / 2;
	double t_deltaTH = maxdist / 8;
	
	t_clusterNum = 0;
	cl.resize(data.size(), -1);
	icenter_of_class.clear();
	for(int i = 0; i < data.size(); i ++)
	{
		if (t_rho[i] > t_rhoTH && t_delta[i] > t_deltaTH)
		{
			cl[i] = t_clusterNum;
			icenter_of_class.push_back(i);
			t_clusterNum++;
		}
	}

	cout << "Performing assignment of each points." << endl;
	for(int i = 0; i < data.size(); i ++)
	{
		if(cl[t_orderrho[i].id] == -1)
		{
			cl[t_orderrho[i].id] = cl[t_neighbor[t_orderrho[i].id]];
		}
	}


	for(int i = 0;i < data.size(); i ++)
	{
		data[i].label = cl[i];
	}
	for(int i = 0;i < icenter_of_class.size(); i ++)
	{
		int centerIDofclsI = icenter_of_class[i];
		data[centerIDofclsI].center = true;
	}
}

void Cluster::Output(string filename)
{
    ofstream fout(filename);
    if(! fout)
    {
    	cout << "WRITING FILE ERROR" << endl;
    	return;
	}
	for(int i = 0; i < data.size(); i ++) fout << data[i].label << endl;
	
	fout.close();
}

int main()
{
	string Input_file = "datasets\\Spiral_cluster=3.txt";
	string Output_file = "result\\DenPeak_Spiral_cluster=3.txt";
	Cluster c(Input_file);
	c.getdist();
	c.getdc();
	c.getrho();
	c.getdelta();
	c.assign();
	c.Output(Output_file);		
}

