#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "maths/stats.h"

//#define UP_ARROW "\033[31m↑\033[0m"
//#define DOWN_ARROW "\033[34m↓\033[0m"

#define UP_ARROW "0"
#define DOWN_ARROW "1"

class RandomService
{
private:
    random_device rd;
protected:
    mt19937 eng;
public:
    RandomService() : eng(rd()) {}
    
    double throwarandompoint_normal(double mean=0,double stdev=1.0) //PCMODIFIED: ADDITION (added mean and stdev parameters)
    {
        normal_distribution<double> distr(mean,stdev);
        return distr(eng);
    }
    double throwarandompoint(double min=0, double max=1)
    {
        uniform_real_distribution<> distr(min,max);
        return distr(eng);
    }
	
	inline bool toss(double wt=0.5) {return (throwarandompoint(0,1)<=wt);}
	inline unsigned long randInt(unsigned long int L) {return throwarandompoint(0,L);}
	template<class T> const T& choice(const std::vector<T>& v) {return v[(int)throwarandompoint(0,v.size())];}
	template<class T> T& choice(std::vector<T>& v) {return v[(int)throwarandompoint(0,v.size())];}
	int weightSelect(double* A,int n)
	{
		double maxV=sum(A,n);
		double rn=throwarandompoint(0,1)*maxV;
		int i;
		for(i=0;i<n;i++)
		{
			rn-=A[i];
			if(rn<0) break;
		}
		return i;
	}
	int weightSelect(const std::vector<double>& wts)
	{
		double v[wts.size()]; for(int i=0;i<wts.size();i++) v[i]=wts[i];
		return weightSelect(v,wts.size());
	}
};

template<typename O, int N> class OverlaidIsingGrid
{
protected:
    bool** spins;
    O** overlay;
    RandomService* rserv=nullptr;
    bool serves_self=true;
public:
    OverlaidIsingGrid()
    {
        spins=new bool*[N];
        overlay=new O*[N];
        for(int i=0;i<N;i++)
        {
            spins[i]=new bool[N];
            overlay[i]=new O[N];
        }
        this->initializeRandom();
    }
    ~OverlaidIsingGrid()
    {
        for(int i=0;i<N;i++)
        {
            delete[] spins[i];
            delete[] overlay[i];
        }
        delete[] spins;
        delete[] overlay;
        if(serves_self) delete rserv;
    }
    
    inline void bindRandomService(RandomService& r)
    {
        if(serves_self) delete rserv;
        rserv=&r;
        serves_self=false;
    }
    
    inline void initializeRandomService()
    {
        if(serves_self) delete rserv;
        rserv=new RandomService();
    }
    
    void initializeRandom()
    {
        if(!rserv) this->initializeRandomService();
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++) spins[i][j]=rserv->toss();
        }
    }
    void resetOverlay(const O& v)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++) overlay[i][j]=v;
        }
    }
    
    bool** cloneSpinMatrix()const
    {
        bool** ret=new bool*[N];
        for(int i=0;i<N;i++) ret[i]=new bool[N];
        
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++) ret[i][j]=this->spins[i][j];
        }
        return ret;
    }
    
    void writeToStream(bool with_mirrors=false, std::ostream& os=std::cout)const
    {
        const std::string SPACE_STRING=" ";
        if(with_mirrors)
        {
            os << " "+SPACE_STRING+"|"+SPACE_STRING;
            for(int j=0;j<N;j++) os << ((this->spins[N-1][j])?UP_ARROW:DOWN_ARROW) <<SPACE_STRING;
            os << "|\n";
            os << "--- "; for(int j=0;j<N;j++) os <<"-"<<SPACE_STRING; os <<"---\n";
        }
        for(int i=0;i<N;i++)
        {
            if(with_mirrors) os << (this->spins[i][N-1]?UP_ARROW:DOWN_ARROW) << SPACE_STRING+"|"+SPACE_STRING;
            for(int j=0;j<N;j++)
                os << ((this->spins[i][j])?UP_ARROW:DOWN_ARROW) <<SPACE_STRING;
            if(with_mirrors) os <<"|"+SPACE_STRING<< (this->spins[i][0]?UP_ARROW:DOWN_ARROW)<<"";
            os << "\n";
        }
        if(with_mirrors)
        {
            os << "--- "; for(int j=0;j<N;j++) os <<"-"<<SPACE_STRING; os <<"---\n";
            os << " "+SPACE_STRING+"|"+SPACE_STRING;
            for(int j=0;j<N;j++) os << ((this->spins[0][j])?UP_ARROW:DOWN_ARROW) <<SPACE_STRING;
            os << "|\n";
        }
    }
};

template<int N> class ClusterFlipIsingModel : public OverlaidIsingGrid<int,N>
{
protected:
    double Jx=1.0,Jy=1.0;
    double beta=1.0;
    int max_cluster=-1;
    std::vector<bool**> bondingGrids; // 4 grids (up, right, down, left - in that order)
public:
    ClusterFlipIsingModel(double B=1.0, double jx=1.0, double jy=1.0) : OverlaidIsingGrid<int,N>() {Jx=jx; Jy=jy; beta=B; _initializeBondingGrids();}
    ~ClusterFlipIsingModel() {_deleteBondingGrids();}
    
protected:
    void _initializeBondingGrids()
    {
        for(int I=0;I<4;I++)
        {
            bondingGrids.push_back(new bool*[N]);
            for(int i=0;i<N;i++)
            {
                bondingGrids[I][i]=new bool[N];
                for(int j=0;j<N;j++) bondingGrids[I][i][j]=false;
            }
        }
    }
    void _deleteBondingGrids()
    {
        for(int I=0;I<4;I++)
        {
            for(int i=0;i<N;i++) delete[] bondingGrids[I][i];
            delete[] bondingGrids[I];
        }
    }
public:
    inline void resetOverlay() {OverlaidIsingGrid<int,N>::resetOverlay(-1);}
    inline bool xthrow()const {return this->rserv->toss(1.0-::exp(-2.0*beta*Jx));}
    inline bool ythrow()const {return this->rserv->toss(1.0-::exp(-2.0*beta*Jy));}
    
    void constructBondingGrids()
    {
        if(!this->rserv) this->initializeRandomService();
        for(int i=0;i<N;i++)
        {
            int iprev=(i-1); if(iprev<0) iprev=N-1; //,inext=(i+1)%N;
            for(int j=0;j<N;j++)
            {
                int jprev=(j-1); if(jprev<0) jprev=N-1; // jnext=(j+1)%N;
                
                if(this->spins[i][j]==this->spins[iprev][j] && this->ythrow())
                {
                    bondingGrids[0][i][j]=true;
                    bondingGrids[2][iprev][j]=true;
                }
                else
                {
                    bondingGrids[0][i][j]=false;
                    bondingGrids[2][iprev][j]=false;
                }
                
                if(this->spins[i][j]==this->spins[i][jprev] && this->xthrow())
                {
                    bondingGrids[3][i][j]=true;
                    bondingGrids[1][i][jprev]=true;
                }
                else
                {
                    bondingGrids[3][i][j]=false;
                    bondingGrids[1][i][jprev]=false;
                }
            }
        }
    }
    
    void _updateClusterLabels(int i,int j, int rid)
    {
        this->overlay[i][j]=rid;
        if(this->bondingGrids[0][i][j])
        {
            int iprev=i-1; if(iprev<0) iprev=N-1;
            if(this->overlay[iprev][j]<0) this->_updateClusterLabels(iprev,j,rid);
        }
        if(this->bondingGrids[3][i][j])
        {
            int jprev=j-1; if(jprev<0) jprev=N-1;
            if(this->overlay[i][jprev]<0) this->_updateClusterLabels(i,jprev,rid);
        }
        if(this->bondingGrids[2][i][j])
        {
            int iprev=i+1; if(iprev>=N) iprev=0;
            if(this->overlay[iprev][j]<0) this->_updateClusterLabels(iprev,j,rid);
        }
        if(this->bondingGrids[1][i][j])
        {
            int jprev=j+1; if(jprev>=N) jprev=0;
            if(this->overlay[i][jprev]<0) this->_updateClusterLabels(i,jprev,rid);
        }
    }
    
    void updateClusterLabels()
    {
        this->resetOverlay();
        int running_id=0;
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                if(this->overlay[i][j]<0) this->_updateClusterLabels(i,j,running_id++);
            }
        }
        max_cluster=running_id;
    }
    
    inline std::string bondingSymbol(int i,int j)const
    {
        const std::vector<std::string> contact_chars={"o",  //0
            u8"↑", u8"→", u8"└", u8"↓",  // 1, 2, (1+2), 4
            u8"|", u8"┌", u8"├", u8"←", // (1+4), (2+4), (1+2+4), 8
            u8"┘", u8"–", u8"┴", u8"┐", // (8+1), (8+2), (1+2+8), (4+8)
            u8"┤", u8"┬", u8"+" // (8+4+1), (8+4+2), (8+4+2+1)
        };
        short int contact_code=0;
        
        if(this->bondingGrids[0][i][j]) contact_code+=1;
        if(this->bondingGrids[1][i][j]) contact_code+=2;
        if(this->bondingGrids[2][i][j]) contact_code+=4;
        if(this->bondingGrids[3][i][j]) contact_code+=8;
        return contact_chars[contact_code];
    }
    void bondingGridToStream(std::ostream& os=std::cout)const
    {
        const std::string SPACE_STRING=" ";
        for(int i=0;i<N;i++)
        {
            //if(with_mirrors) os << (this->spins[i][N-1]?UP_ARROW:DOWN_ARROW) << SPACE_STRING+"|"+SPACE_STRING;
            for(int j=0;j<N;j++)
                os << this->bondingSymbol(i,j) <<SPACE_STRING;
            //if(with_mirrors) os <<"|"+SPACE_STRING<< (this->spins[i][0]?UP_ARROW:DOWN_ARROW)<<"";
            os << "\n";
        }
    }
    
    void clusterLabelsToStream(std::ostream& os=std::cout)const
    {
        const std::string SPACE_STRING=" ";
        char* buf=new char[8];
        std::string fmt="%02d";
        for(int i=0;i<N;i++)
        {
            //if(with_mirrors) os << (this->spins[i][N-1]?UP_ARROW:DOWN_ARROW) << SPACE_STRING+"|"+SPACE_STRING;
            for(int j=0;j<N;j++)
            {
                sprintf(buf,fmt.c_str(),this->overlay[i][j]);
                os << buf <<SPACE_STRING;
            }
            //if(with_mirrors) os <<"|"+SPACE_STRING<< (this->spins[i][0]?UP_ARROW:DOWN_ARROW)<<"";
            os << "\n";
        }
        delete[] buf;
    }
    
    void step()
    {
        this->constructBondingGrids();
        this->updateClusterLabels();
        int selclust=(int)(this->rserv->throwarandompoint(0,1)*max_cluster-1e-6);
        //std::cout << "Selected cluster: "<<selclust<<" to flip\n";
        if(this->rserv->toss())
        {
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    if(this->overlay[i][j]==selclust) this->spins[i][j]=!this->spins[i][j];
                }
            }
        }
    }
};
