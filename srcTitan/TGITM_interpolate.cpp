//$Id$
#include <iostream>
#include <fstream>
#include <sstream>
//#include <array>
#include "Titan.h"


using namespace std;

	double Titan::tgitm_exobase::tgitm_grid[Titan::tgitm_exobase::ndes][10], Titan::tgitm_exobase::interp_val[7];
	
	
	void Titan::tgitm_exobase::read_tgitm() {
	const int nLatitudes = 72;

	
	cout<<"!!!!!!!!!!!!!"<<"\t"<<Titan::tgitm_exobase::ndes<<endl;
	ifstream fin("TGITM_EXOBASE.dat");
	
	if(fin.is_open())
	{
		for(int i=0;i<Titan::tgitm_exobase::ndes;i++)
		{
			for(int j=0;j<10;j++)
			{
				fin >> Titan::tgitm_exobase::tgitm_grid[i][j];
				//cout<<"!!!!!!!!!!!!!"<<"\t"<<Titan::tgitm_exobase::tgitm_grid[i][j]<<endl;
			}
		}
	}
	else
	{
			cout<<"THE FILE DID NOT OPEN"<<endl;
	}
	
	}


  void Titan::tgitm_exobase::tgitm_interpolate(double polar, double azimuth)
//Angles should be in Titan center frame
{
	double lat, lng, fp[4];
	int i_lat, i_lng, i_cell, itp_point[4];
	const double tgitm_lb_lat=1.25, tgitm_hb_lat=178.75, 
	tgitm_lb_lng=2.5, tgitm_hb_lng=357.5, dlat=2.5, dlng=5.0;
	const int nLatitudes = 72;
	const double rad2dg=57.2957795;

	//lat=polar;
	//lng=azimuth;
	
	//Spherical coordinates conveted to lattitude and east longitude for interpolation
	lat=180.0-polar*rad2dg;
	lng=360.0-azimuth*rad2dg;
	
	if(lng < tgitm_lb_lng) lng=tgitm_lb_lng+dlng*0.25;
	if(lng > tgitm_hb_lng) lng=tgitm_hb_lng-dlng*0.25;
	if(lat > tgitm_hb_lat) lat=tgitm_hb_lat-dlat*0.25;
	if(lat < tgitm_lb_lat) lat=tgitm_lb_lat+dlat*0.25;
	
	i_lat=int((lat-tgitm_lb_lat)/dlat);
	i_lng=int((lng-tgitm_lb_lng)/dlng);
	i_cell=i_lng*(nLatitudes-1)+i_lat;

	itp_point[0]=i_cell+i_lng;
	itp_point[1]=itp_point[0]+1;
	itp_point[2]=itp_point[0]+nLatitudes;
	itp_point[3]=itp_point[2]+1;
	
	
	for(int i=0; i<7; i++)
	{
		fp[0] = Titan::tgitm_exobase::tgitm_grid[itp_point[0]][i+2];
		fp[1] = Titan::tgitm_exobase::tgitm_grid[itp_point[1]][i+2];
		fp[2] = Titan::tgitm_exobase::tgitm_grid[itp_point[2]][i+2];
		fp[3] = Titan::tgitm_exobase::tgitm_grid[itp_point[3]][i+2];
	
		Titan::tgitm_exobase::interp_val[i]=fp[0]*(Titan::tgitm_exobase::tgitm_grid[itp_point[2]][0]-lng)*(Titan::tgitm_exobase::tgitm_grid[itp_point[3]][1]-lat)+
		
		fp[1]*(Titan::tgitm_exobase::tgitm_grid[itp_point[2]][0]-lng)*(lat -Titan::tgitm_exobase::tgitm_grid[itp_point[2]][1])+
		
		fp[2]*(lng-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][0])*(Titan::tgitm_exobase::tgitm_grid[itp_point[1]][1]-lat)+
		
		fp[3]*(lng-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][0])*(lat-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][1]);
			
		Titan::tgitm_exobase::interp_val[i]=Titan::tgitm_exobase::interp_val[i]/(Titan::tgitm_exobase::tgitm_grid[itp_point[2]][0]-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][0])
		/(Titan::tgitm_exobase::tgitm_grid[itp_point[1]][1]-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][1]);
		
		
	}
	/*cout<<"start"<<endl;
	cout<<itp_point[0]<<"\t"<<itp_point[1]<<"\t"<<itp_point[2]<<"\t"<<itp_point[3]<<"\t"<<endl;
	cout<<fp[0]<<"\t"<<fp[1]<<"\t"<<fp[2]<<"\t"<<fp[3]<<"\t"<<endl;
	cout<<"lat"<<"\t"<<lat<<"\t"<<"lng"<<"\t"<<lng<<"\t"<<Titan::tgitm_exobase::interp_val[3]<<"\t"<<endl;
	cout<<"end"<<endl<<endl<<endl;*/
}


/*int main()
{
	const int ndes=5184;
	double tgitm_grid[ndes][7], interp_val[4];
	double theta,phi,dtheta,dphi;
	int ntheta=100,nphi=100,ii[4],iii,cnt;
	cout<<"Hello"<<endl;
	
	read_tgitm(tgitm_grid, ndes);
	
	for(int i=0;i<ndes;i++)
		{
			for(int j=0;j<6;j++)
			{
				cout << tgitm_grid[i][j]<<'\t';
			}
			cout <<endl;
		}
		dtheta=360.0/ntheta;
		dphi=180.0/nphi;
		theta=0.0;
		phi=0.0;
		ofstream fo("TGITM_cc.dat");
		fo<< "VARIABLES= "<<"long, "<<"lat, "<<"n[N2], " <<"n[CH4], "<<"n{H2], "<<"T[K] "<<endl;
		fo<<"ZONE NODES= "<<ntheta*nphi<<","<<"ELEMENTS= "<< (ntheta-1)*(nphi-1)<<", "<< "DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<<endl;
		
		
		
		for(int i=0;i<ntheta;i++)
		{
			for( int j=0;j<nphi;j++)
			{
				tgitm_interpolate(phi, theta, interp_val, tgitm_grid);
				fo<<theta<<"\t"<<phi<<"\t"<<interp_val[0]<<"\t"<<interp_val[1]<<"\t"<<
				interp_val[2]<<"\t"<<interp_val[3]<<"\t"<<endl;
				phi=phi+dphi;
				//return 0;
			}
			phi=0.0;
			theta=theta+dtheta;
		}
		iii=0;
		cnt=0;
		for(int i=1;i<=(ntheta-1)*(nphi-1);i++)
		{
			cnt=cnt+1;
			if(cnt % nphi != 0)
			{
				iii=iii+1;
			}
			else
			{
				cnt=1;
				iii=iii+2;
			}
			
			ii[0]=iii;
			ii[1]=ii[0]+1;
			ii[2]=ii[1]+nphi;
			ii[3]=ii[2]-1;
			fo<<ii[0]<<"\t"<<ii[1]<<"\t"<<ii[2]<<"\t"<<ii[3]<<"\t"<<endl;
			
		}
		fo.close();
		
		
		return 0;
		
}
*/
