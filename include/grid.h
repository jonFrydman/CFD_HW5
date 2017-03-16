#ifndef GRID_H
#define GRID_H
#include<fstream>
#include<string>
#include<cstdlib>
#include<algorithm>
#include<vector>
#include<cmath>
#include<iostream>
using namespace std;
class grid{
//Reads external grid files, computes additional grid data, holds geometric data
//Any particular flow problem should only need one grid
    typedef std::vector< std::vector<double> > vec2D;
    public:
        grid(std::string xfilename, std::string yfilename){
            loadGrid(xfilename, yfilename); //sets x and y corner points
            defineCenterPoints();
            defineSideLengths();
            defineSideNormals();
            defineAreas();
            defineIJNormals();
			defineDeltaS();
        }
        vec2D xCorner, yCorner, xCenter, yCenter;
        vec2D xSside, ySside, xWside, yWside;
        vec2D xSnorm, ySnorm, xWnorm, yWnorm;
        vec2D xSdeltas, ySdeltas, xWdeltas, yWdeltas;
        vec2D xInorm, yInorm, xJnorm, yJnorm; //the centered, averaged normals
        vec2D area;
        int N, M; //number of points in the xi, and eta directions
    private:
        double **dx, **dy; //length parameters of cells between corners. Indices match edge-point indices.
        void loadGrid(std::string xfilename, std::string yfilename){ //Define N and M in here; x_corners = new double[][]; more initializing stuff
            std::ifstream xfile, yfile;
            std::string readin;
            int lineNum=0;
            xfile.open(xfilename.c_str());
            yfile.open(yfilename.c_str());
            if(!xfile.good()){std::cout<<"Error Reading x-grid file!\n";}
            if(!yfile.good()){std::cout<<"Error Reading y-grid file!\n";}

            //Figure out size of the csv table
            getline(xfile, readin);
            lineNum++;
            M=std::count(readin.begin(),readin.end(),',')+1;
            while(xfile.good()){
                getline(xfile, readin);
                lineNum++;
            }
            N=lineNum-1;
            xfile.clear();
            xfile.seekg(0, xfile.beg);//restart xfile from the beginning

            //Start Loading in data now that we know how much information there is
            //Allocate the correct amount of memory for xCorners and yCorners
            xCorner.resize(N, std::vector<double>(M));
            yCorner.resize(N, std::vector<double>(M));

            for(int i=0; i<N; i++){
                //read all data entries except the last in the line (the last line is the only one that doesnt have a comma at the end of it)

				// for fucks sake. Some of the grid points were negative zeros!!!! 3AM frustration brought to you by Mike Wennerstrom
                for(int j=0; j<M-1; j++){
                    getline(xfile,readin,',');
                    xCorner[i][j]=atof(readin.c_str());
                    getline(yfile,readin,',');
                    yCorner[i][j]= atof(readin.c_str());
                }
                getline(xfile,readin);
                xCorner[i][M-1]=atof(readin.c_str());
                getline(yfile,readin);
                yCorner[i][M-1]= atof(readin.c_str());
            }
            xfile.close();
            yfile.close();

        }
        void defineCenterPoints(){ //interpolate cell centers
            xCenter.resize(N-1, std::vector<double>(M-1));
            yCenter.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    xCenter[i][j]=(xCorner[i][j]+xCorner[i][j+1]+xCorner[i+1][j]+xCorner[i+1][j+1])/4;
                    yCenter[i][j]=(yCorner[i][j]+yCorner[i][j+1]+yCorner[i+1][j]+yCorner[i+1][j+1])/4;
                }
            }
        }
        void defineAreas(){
            area.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    area[i][j]=0.5*abs(xSside[i][j]*yWside[i][j]-xWside[i][j]*ySside[i][j]);
                }
            }
        }
        void defineSideLengths(){
            xSside.resize(N-1, std::vector<double>(M)); // doesn't bother past N-1
            ySside.resize(N-1, std::vector<double>(M));
            xWside.resize(N, std::vector<double>(M-1)); // doesn't have to look past M-1
            yWside.resize(N, std::vector<double>(M-1));
            //Side length on the south side and west side of cell
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    xSside[i][j]=abs(xCorner[i+1][j]-xCorner[i][j]);
                    ySside[i][j]=abs(yCorner[i+1][j]-yCorner[i][j]);
                    xWside[i][j]=abs(xCorner[i][j+1]-xCorner[i][j]);
                    yWside[i][j]=abs(yCorner[i][j+1]-yCorner[i][j]);

					xWside[N - 1][j] = abs(xCorner[0][j + 1] - xCorner[0][j]); // boundary conditions
					yWside[N - 1][j] = abs(yCorner[0][j + 1] - yCorner[0][j]);

                }
                xSside[i][M-1]=abs(xCorner[i+1][M-1]-xCorner[i][M-1]);
                ySside[i][M-1]=abs(yCorner[i+1][M-1]-yCorner[i][M-1]); // checked at 3/15 9pm
            }
        }
        void defineSideNormals(){
            xSnorm.resize(N-1, std::vector<double>(M));
            ySnorm.resize(N-1, std::vector<double>(M));
            xWnorm.resize(N-1, std::vector<double>(M-1));
            yWnorm.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    double Slength=sqrt(pow(xSside[i][j],2)+pow(ySside[i][j],2));
                    double Wlength=sqrt(pow(xWside[i][j],2)+pow(yWside[i][j],2));
                    ySnorm[i][j]=(xCorner[i+1][j]-xCorner[i][j])/Slength; //checked 3/15 at 9pm
                    xSnorm[i][j]=-(yCorner[i+1][j]-yCorner[i][j])/Slength;
                    yWnorm[i][j]=-(xCorner[i][j+1]-xCorner[i][j])/Wlength;
                    xWnorm[i][j]=(yCorner[i][j+1]-yCorner[i][j])/Wlength;
                }
                double Slength=sqrt(pow(xSside[i][M-1],2)+pow(ySside[i][M-1],2));
                ySnorm[i][M-1]=(xCorner[i+1][M-1]-xCorner[i][M-1])/Slength;
                xSnorm[i][M-1]=-(yCorner[i+1][M-1]-yCorner[i][M-1])/Slength; //checked 3/15 at 11:45
            }
        }
        void defineDeltaS(){
            xSdeltas.resize(N-1, std::vector<double>(M));
            ySdeltas.resize(N-1, std::vector<double>(M));
            xWdeltas.resize(N-1, std::vector<double>(M-1));
            yWdeltas.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    ySdeltas[i][j]=(xCorner[i+1][j]-xCorner[i][j]);
                    xSdeltas[i][j]=-(yCorner[i+1][j]-yCorner[i][j]);
                    yWdeltas[i][j]=-(xCorner[i][j+1]-xCorner[i][j]);
                    xWdeltas[i][j]=(yCorner[i][j+1]-yCorner[i][j]);
                }
                double Slength=sqrt(pow(xSside[i][M-1],2)+pow(ySside[i][M-1],2));
                ySdeltas[i][M-1]=(xCorner[i+1][M-1]-xCorner[i][M-1]); // fixed from y=-x, x=y to y=x and x=-y 3/15 9pm
                xSdeltas[i][M-1]=-(yCorner[i+1][M-1]-yCorner[i][M-1]);
            }
        }
        void defineIJNormals(){
            xInorm.resize(N-1, std::vector<double>(M-1));
            yInorm.resize(N-1, std::vector<double>(M-1));
            xJnorm.resize(N-1, std::vector<double>(M-1));
            yJnorm.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-2; i++){ // Altered to include row i = 127 and column j=63 instead of all zeros
                for(int j=0; j<M-2; j++){
                    xInorm[i][j]=0.5*(xWnorm[i][j]+xWnorm[i+1][j]);
                    yInorm[i][j]=0.5*(yWnorm[i][j]+yWnorm[i+1][j]);
                    xJnorm[i][j]=0.5*(xSnorm[i][j]+xSnorm[i][j+1]);
                    yJnorm[i][j]=0.5*(ySnorm[i][j]+ySnorm[i][j+1]);
                }
            }

			//To get the boundary values for these s's along the cells in i = 127. Note that the iplus i + 1 values are replaced with 0
			for (int j = 0; j < M - 2; j++) {

				int i = N - 2;

				double Slength = sqrt(pow(xSside[i][j], 2) + pow(ySside[i][j], 2));
				double Slength_jplus = sqrt(pow(xSside[i][j+1], 2) + pow(ySside[i][j+1], 2));
				double Wlength = sqrt(pow(xWside[i][j], 2) + pow(yWside[i][j], 2));
				double Wlength_iplus = sqrt(pow(xWside[i+1][j], 2) + pow(yWside[i+1][j], 2));

				double temp_ySnorm = (xCorner[i + 1][j] - xCorner[i][j]) / Slength;
				double temp_ySnorm_jplus = (xCorner[i + 1][j+1] - xCorner[i][j+1]) / Slength_jplus;

				double temp_xSnorm = -(yCorner[i + 1][j] - yCorner[i][j]) / Slength;
				double temp_xSnorm_jplus = -(yCorner[i + 1][j+1] - yCorner[i][j+1]) / Slength_jplus;

				double temp_yWnorm = -(xCorner[i][j + 1] - xCorner[i][j]) / Wlength;
				double temp_yWnorm_iplus = -(xCorner[i+1][j + 1] - xCorner[i+1][j]) / Wlength_iplus;

				double temp_xWnorm = (yCorner[i][j + 1] - yCorner[i][j]) / Wlength;
				double temp_xWnorm_iplus = (yCorner[i+1][j + 1] - yCorner[i+1][j]) / Wlength_iplus;

				xInorm[i][j] = 0.5*(temp_xWnorm + temp_xWnorm_iplus);
				yInorm[i][j] = 0.5*(temp_yWnorm + temp_yWnorm_iplus);
				xJnorm[i][j] = 0.5*(temp_xSnorm + temp_xSnorm_jplus);
				yJnorm[i][j] = 0.5*(temp_ySnorm + temp_ySnorm_jplus);

			}

			//To get the boundary values for these s's along the cells in j = 63. Note the boundary condition when i = 127
			for (int i = 0; i < N - 1; i++) {

				int j = M - 2;

				double Slength = sqrt(pow(xSside[i][j], 2) + pow(ySside[i][j], 2));
				double Slength_jplus = sqrt(pow(xSside[i][j + 1], 2) + pow(ySside[i][j + 1], 2));
				double Wlength = sqrt(pow(xWside[i][j], 2) + pow(yWside[i][j], 2));
				double Wlength_iplus = sqrt(pow(xWside[i + 1][j], 2) + pow(yWside[i + 1][j], 2));

				double temp_ySnorm = (xCorner[i + 1][j] - xCorner[i][j]) / Slength;
				double temp_ySnorm_jplus = (xCorner[i + 1][j + 1] - xCorner[i][j + 1]) / Slength_jplus;

				double temp_xSnorm = -(yCorner[i + 1][j] - yCorner[i][j]) / Slength;
				double temp_xSnorm_jplus = -(yCorner[i + 1][j + 1] - yCorner[i][j + 1]) / Slength_jplus;

				double temp_yWnorm = -(xCorner[i][j + 1] - xCorner[i][j]) / Wlength;
				double temp_yWnorm_iplus = -(xCorner[i + 1][j + 1] - xCorner[i + 1][j]) / Wlength_iplus;

				double temp_xWnorm = (yCorner[i][j + 1] - yCorner[i][j]) / Wlength;
				double temp_xWnorm_iplus = (yCorner[i + 1][j + 1] - yCorner[i + 1][j]) / Wlength_iplus;


				xInorm[i][j] = 0.5*(temp_xWnorm + temp_xWnorm_iplus); // 3/15 this needs to be fixed
				yInorm[i][j] = 0.5*(temp_yWnorm + temp_yWnorm_iplus);

				xJnorm[i][j] = 0.5*(temp_xSnorm + temp_xSnorm_jplus);
				yJnorm[i][j] = 0.5*(temp_ySnorm + temp_ySnorm_jplus);

			}

        }
};

#endif // GRID_H
